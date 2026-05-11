#!/usr/bin/env python3

import subprocess, sys, shutil
from pathlib import Path

BASE = Path("/project/home/p201228/gpcr_screening")

# Cadena del receptor para cada PDB
RECEPTOR_CHAIN = {
    "8J6P": "R",
    "8J6Q": "R",
    "8JZ7": "A",
    "9IZC": "A",
    "8IJA": "A",
    "8J6J": "R",
    "8JEF": "A",
    "8JEI": "A",
    "9JID": "A",
}

RECEPTORS = {
    "GPR109A": ["8J6P", "8J6Q", "8JZ7", "9IZC", "8IJA", "8J6J"],
    "GPR109B": ["8JEF", "8JEI", "9JID"],
}

AD4_LIG_TYPES = "A C HD N NA OA SA S Cl F Br I".split()

# Espaciado de la cuadrícula: 0.500 Å para docking ciego
GRID_SPACING   = 0.500

# Margen alrededor del receptor (Å)
PROTEIN_MARGIN = 5.0

# 1. Extraer cadena del receptor del PDB original

def extract_chain(pdb_id, raw_pdb, out_pdb, chain):
    if out_pdb.exists() and out_pdb.stat().st_size > 0:
        print(f"  [SKIP] {pdb_id} cadena {chain} ya extraída"); return True

    n_atoms = 0
    with open(raw_pdb) as fin, open(out_pdb, 'w') as fout:
        for line in fin:
            record = line[:6].strip()
            if record in ("ATOM", "HETATM"):
                # Columna 22 (índice 21) = chain ID en formato PDB estándar
                line_chain = line[21]
                if line_chain == chain:
                    fout.write(line)
                    n_atoms += 1
            elif record in ("TER", "END", "REMARK", "HEADER", "TITLE"):
                fout.write(line)

    if n_atoms == 0:
        print(f"  [ERROR] {pdb_id}: no se encontraron átomos en cadena {chain}",
              file=sys.stderr)
        out_pdb.unlink(missing_ok=True)
        return False

    print(f"  [OK] {pdb_id} cadena {chain} extraída ({n_atoms} átomos)")
    return True

# 2. Limpiar con pdbfixer (solo eliminar agua y añadir H)
def clean_chain(pdb_id, in_f, out_f):
    if out_f.exists() and out_f.stat().st_size > 0:
        print(f"  [SKIP] {pdb_id} limpieza ya hecha"); return True
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
        fixer = PDBFixer(filename=str(in_f))
        fixer.removeHeterogens(keepWater=False)
        fixer.addMissingHydrogens(7.4)
        with open(out_f, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
        print(f"  [OK] {pdb_id} limpiado"); return True
    except Exception as e:
        # Fallback: limpiar manualmente con grep (para HIS problemáticas)
        print(f"  [WARN] {pdb_id} pdbfixer: {e} — usando fallback grep")
        try:
            with open(in_f) as fin, open(out_f, 'w') as fout:
                for line in fin:
                    if line.startswith("ATOM"):
                        fout.write(line)
            print(f"  [OK-FB] {pdb_id} limpiado (solo ATOM)"); return True
        except Exception as e2:
            print(f"  [ERROR] {pdb_id}: {e2}", file=sys.stderr); return False

# 3. Convertir a PDBQT

def to_pdbqt(pdb_id, in_f, out_f):
    if out_f.exists() and out_f.stat().st_size > 0:
        print(f"  [SKIP] {pdb_id} PDBQT ya existe"); return True

    r = subprocess.run(
        ["mk_prepare_receptor.py", "-i", str(in_f), "-o", str(out_f)],
        capture_output=True, text=True
    )
    if r.returncode == 0 and out_f.exists() and out_f.stat().st_size > 0:
        print(f"  [OK] {pdb_id} → PDBQT (meeko)"); return True

    # Fallback obabel
    r2 = subprocess.run(
        ["obabel", str(in_f), "-O", str(out_f),
         "--partialcharge", "gasteiger", "-h"],
        capture_output=True, text=True
    )
    if r2.returncode == 0 and out_f.exists() and out_f.stat().st_size > 0:
        print(f"  [OK-OB] {pdb_id} → PDBQT (obabel)"); return True

    print(f"  [ERROR] {pdb_id} PDBQT fallido\n  {r.stderr[:200]}", file=sys.stderr)
    return False

# 4. Calcular grid box ciego sobre la cadena del receptor

def get_rec_types(pdbqt_file):
    types = set()
    with open(pdbqt_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                if len(parts) >= 12:
                    types.add(parts[-1].strip())
    return sorted(types)

def make_gpf_blind_chain(pdb_id, chain, pdbqt_file, gpf_out,
                          spacing=GRID_SPACING, margin=PROTEIN_MARGIN):

    rec_types = get_rec_types(pdbqt_file)
    print(f"  [INFO] {pdb_id} receptor_types: {' '.join(rec_types)}")

    coords = []
    with open(pdbqt_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    coords.append((
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54])
                    ))
                except ValueError:
                    continue

    if not coords:
        print(f"  [ERROR] {pdb_id}: sin coordenadas en PDBQT", file=sys.stderr)
        return False

    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]

    # Centro geométrico de la cadena del receptor
    cx = (max(xs) + min(xs)) / 2
    cy = (max(ys) + min(ys)) / 2
    cz = (max(zs) + min(zs)) / 2

    # Dimensiones: rango de la cadena + margen
    rx = max(xs) - min(xs) + 2 * margin
    ry = max(ys) - min(ys) + 2 * margin
    rz = max(zs) - min(zs) + 2 * margin

    def to_even(size, sp=spacing):
        pts = int(size / sp)
        return pts if pts % 2 == 0 else pts + 1

    nx, ny, nz = to_even(rx), to_even(ry), to_even(rz)

    print(f"  [BLIND-CHAIN] {pdb_id} (cadena {chain}):")
    print(f"    Centro receptor : ({cx:.2f}, {cy:.2f}, {cz:.2f}) Å")
    print(f"    Tamaño caja     : {nx*spacing:.1f} × {ny*spacing:.1f} × {nz*spacing:.1f} Å")
    print(f"    Puntos grid     : {nx} × {ny} × {nz} = {nx*ny*nz:,} puntos")
    print(f"    Espaciado       : {spacing} Å  |  Margen: {margin} Å")
    print(f"    Átomos receptor : {len(coords)}")

    if max(nx, ny, nz) > 200:
        print(f"  [WARN] Caja >200 pts en alguna dimensión — considera spacing=0.6")

    fld_name = f"{pdb_id}_blind_chain{chain}.maps.fld"
    with open(gpf_out, 'w') as f:
        f.write(f"# GPF — docking ciego cadena {chain} — {pdb_id}\n")
        f.write(f"# Solo receptor GPR (cadena {chain}, {len(coords)} átomos)\n")
        f.write(f"npts {nx} {ny} {nz}\n")
        f.write(f"gridfld {fld_name}\n")
        f.write(f"spacing {spacing}\n")
        f.write(f"receptor_types {' '.join(rec_types)}\n")
        f.write(f"ligand_types {' '.join(AD4_LIG_TYPES)}\n")
        f.write(f"receptor {pdbqt_file.name}\n")
        f.write(f"gridcenter {cx:.4f} {cy:.4f} {cz:.4f}\n")
        f.write(f"smooth 0.5\n")
        for at in AD4_LIG_TYPES:
            f.write(f"map {pdb_id}_blind_chain{chain}.{at}.map\n")
        f.write(f"elecmap {pdb_id}_blind_chain{chain}.e.map\n")
        f.write(f"dsolvmap {pdb_id}_blind_chain{chain}.d.map\n")
        f.write(f"dielectric -0.1465\n")
    return True

# 5. Calcular maps con autogrid4

def run_autogrid(pdb_id, chain, gpf, maps_dir):
    fld = maps_dir / f"{pdb_id}_blind_chain{chain}.maps.fld"
    if fld.exists() and fld.stat().st_size > 0:
        print(f"  [SKIP] {pdb_id} maps ya calculados"); return True

    print(f"  Calculando grid maps (puede tardar varios minutos)...")
    r = subprocess.run(
        ["autogrid4", "-p", gpf.name,
         "-l", f"{pdb_id}_blind_chain{chain}.glg"],
        cwd=maps_dir, capture_output=True, text=True
    )
    n = len(list(maps_dir.glob(f"{pdb_id}_blind_chain{chain}.*.map")))
    if fld.exists() and n > 0:
        print(f"  [OK] {pdb_id}: {n} maps generados"); return True

    print(f"  [ERROR] {pdb_id} autogrid4 falló", file=sys.stderr)
    glg = maps_dir / f"{pdb_id}_blind_chain{chain}.glg"
    if glg.exists():
        print('\n'.join(glg.read_text().splitlines()[-12:]), file=sys.stderr)
    return False

# MAIN

print("=" * 62)
print(" Preparación de receptores — DOCKING CIEGO SELECTIVO POR CADENA")
print(f" Espaciado: {GRID_SPACING} Å  |  Margen: {PROTEIN_MARGIN} Å")
print("=" * 62)
print()
print(" Cadenas del receptor identificadas:")
for pdb, chain in RECEPTOR_CHAIN.items():
    protein = "GPR109A" if pdb in RECEPTORS["GPR109A"] else "GPR109B"
    print(f"   {pdb} ({protein}) → cadena {chain}")
print()

ok = fail = 0

for protein, pdbs in RECEPTORS.items():
    print(f"\n════════ {protein} ════════")
    raw_dir  = BASE / protein / "receptors/raw"
    # Directorio específico para docking ciego por cadena
    prep_dir = BASE / protein / "receptors/prepared_chain"
    maps_dir = BASE / protein / "receptors/maps_blind_chain"
    prep_dir.mkdir(parents=True, exist_ok=True)
    maps_dir.mkdir(parents=True, exist_ok=True)

    for pdb_id in pdbs:
        chain = RECEPTOR_CHAIN[pdb_id]
        print(f"\n  ── {pdb_id} (cadena {chain}) ──")

        raw            = raw_dir  / f"{pdb_id}.pdb"
        chain_pdb      = prep_dir / f"{pdb_id}_chain{chain}.pdb"
        clean_pdb      = prep_dir / f"{pdb_id}_chain{chain}_clean.pdb"
        pdbqt          = prep_dir / f"{pdb_id}_chain{chain}.pdbqt"
        maps_pdbqt     = maps_dir / f"{pdb_id}_chain{chain}.pdbqt"
        gpf            = maps_dir / f"{pdb_id}_blind_chain{chain}.gpf"

        if not raw.exists() or raw.stat().st_size == 0:
            print(f"  [FALTA] {raw}"); fail += 1; continue

        # Paso 1: extraer solo la cadena del receptor
        if not extract_chain(pdb_id, raw, chain_pdb, chain):
            fail += 1; continue

        # Paso 2: limpiar (eliminar agua/HETATM, añadir H)
        if not clean_chain(pdb_id, chain_pdb, clean_pdb):
            fail += 1; continue

        # Paso 3: convertir a PDBQT
        if not to_pdbqt(pdb_id, clean_pdb, pdbqt):
            fail += 1; continue

        # Copiar PDBQT al directorio de maps (autogrid lo necesita ahí)
        if not maps_pdbqt.exists():
            shutil.copy(pdbqt, maps_pdbqt)

        # Limpiar maps anteriores si existen
        for old in maps_dir.glob(f"{pdb_id}_blind_chain{chain}.*.map"):
            old.unlink()
        (maps_dir / f"{pdb_id}_blind_chain{chain}.maps.fld").unlink(missing_ok=True)
        gpf.unlink(missing_ok=True)

        # Paso 4: generar GPF con caja sobre la cadena del receptor
        if not make_gpf_blind_chain(pdb_id, chain, maps_pdbqt, gpf):
            fail += 1; continue

        # Paso 5: calcular grid maps
        if not run_autogrid(pdb_id, chain, gpf, maps_dir):
            fail += 1; continue

        ok += 1
        print(f"{pdb_id} completo — solo cadena {chain}")
