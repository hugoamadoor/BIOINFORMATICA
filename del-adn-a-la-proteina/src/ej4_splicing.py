"""Ejercicio 4. Splicing alternativo.

Uso:  python src/ej4_splicing.py                 # gen de juguete con 5 exones
      python src/ej4_splicing.py --ensembl FGFR2  # isoformas reales (requiere Internet)

Parte 1: los exones 1 y 5 son constitutivos y los exones 2, 3 y 4 pueden incluirse
o saltarse, lo que da 2^3 = 8 isoformas. Para cada una se comprueba si conserva la
pauta de lectura, dónde aparece el primer codón de paro y si el ARNm sería
candidato a degradación por NMD (regla de los 50-55 nt).
Parte 2: consulta la API REST de Ensembl y compara los transcritos de un gen humano.
"""
import argparse
import itertools
import json
import time
import urllib.error
import urllib.request

from Bio.Seq import Seq

# Exones del gen de juguete (hebra codificante). Longitudes: 24, 20, 22, 21, 24 nt.
# 20 ≡ 2 y 22 ≡ 1 (mod 3): saltar solo uno de ellos desplaza la pauta de lectura.
# 21 ≡ 0 (mod 3): saltar el exón 4 elimina 7 aminoácidos sin cambiar la pauta.
EXONES = {
    1: "ATGGCCACGGACTTTGCTCAGCGG",
    2: "GATATCGCCCCACGGAATGG",
    3: "AATACATCCTTGTTTACGAATT",
    4: "GTGCGCTATACCATTGGTGGC",
    5: "AGACTGCTTACTCATGTAGTATAA",
}
PROPUESTAS = [(1, 2, 3, 4, 5), (1, 2, 4, 5), (1, 3, 5), (1, 2, 3, 5)]


def analizar(combinacion):
    arnm = "".join(EXONES[e] for e in combinacion)
    prot = str(Seq(arnm[: len(arnm) // 3 * 3]).translate())
    paro = prot.find("*")
    ultima_union = len(arnm) - len(EXONES[combinacion[-1]])
    completa = paro == len(prot) - 1 and len(arnm) % 3 == 0
    info = {
        "exones": "-".join(map(str, combinacion)),
        "nt": len(arnm),
        "en_fase": len(arnm) % 3 == 0,
        "proteina": prot[:paro] if paro >= 0 else prot,
        "completa": completa,
        "sin_paro": paro == -1,
        "paro_prematuro": paro >= 0 and not completa,
    }
    if info["paro_prematuro"]:
        info["dist_union"] = ultima_union - paro * 3   # >0: el paro está antes de la última unión
        info["nmd"] = info["dist_union"] > 50
    return info


def gen_juguete():
    print("Gen de juguete: exón 1 y exón 5 constitutivos; 2, 3 y 4 alternativos.")
    print(f"{'Isoforma':<11}{'nt':>4}  {'Pauta':<9}{'aa':>4}  Proteína / observación")
    todas = [(1, *c, 5) for r in range(4) for c in itertools.combinations((2, 3, 4), r)]
    for comb in sorted(todas, key=lambda c: (-len(c), c)):
        i = analizar(comb)
        marca = " <- propuesta" if comb in PROPUESTAS else ""
        if i["completa"]:
            obs = "proteína completa"
        elif i["sin_paro"]:
            obs = "sin paro en los exones: se leería la 3'UTR (proteína alargada o ARNm non-stop)"
        elif i["nmd"]:
            obs = "paro prematuro >50 nt antes de la última unión: NMD probable"
        elif i["dist_union"] <= 0:
            obs = "paro prematuro en el último exón: escapa al NMD, extremo C alterado"
        else:
            obs = f"paro prematuro solo {i['dist_union']} nt antes de la última unión: escapa al NMD"
        print(f"{i['exones']:<11}{i['nt']:>4}  {'correcta' if i['en_fase'] else 'desplaz.':<9}"
              f"{len(i['proteina']):>4}  {i['proteina']}  ({obs}){marca}")


API = "https://rest.ensembl.org"


def pedir_json(ruta, reintentos=4):
    """GET a la API REST de Ensembl con reintentos ante errores temporales (5xx, 429, red)."""
    url = f"{API}{ruta}{'&' if '?' in ruta else '?'}content-type=application/json"
    for intento in range(1, reintentos + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503, 504) or intento == reintentos:
                raise
            espera = float(e.headers.get("Retry-After", 2 * intento))
        except urllib.error.URLError:
            if intento == reintentos:
                raise
            espera = 2 * intento
        print(f"  (Ensembl no responde, reintento {intento}/{reintentos - 1} en {espera:.0f} s)")
        time.sleep(espera)


def ensembl(simbolo):
    """Lista los transcritos codificantes de un gen humano con su nº de exones y longitud.

    No se usa lookup/symbol?expand=1 sobre el gen completo porque en genes con muchos
    transcritos (FGFR2 tiene 60) la API devuelve a menudo un error 500. En su lugar:
    1) lookup del gen, 2) overlap para listar sus transcritos, 3) lookup de cada uno.
    """
    print(f"Consultando Ensembl ({API}) para {simbolo}...")
    gen = pedir_json(f"/lookup/symbol/homo_sapiens/{simbolo}")
    transcritos = [t for t in pedir_json(f"/overlap/id/{gen['id']}?feature=transcript")
                   if t.get("Parent") == gen["id"]]
    codif = [t for t in transcritos if t.get("biotype") == "protein_coding"]
    print(f"{gen.get('display_name')} ({gen['id']}), cromosoma {gen.get('seq_region_name')}, "
          f"hebra {'+' if gen.get('strand') == 1 else '-'}: {len(transcritos)} transcritos, "
          f"{len(codif)} codificantes de proteína. Descargando detalles...")
    filas = []
    for t in codif:
        detalle = pedir_json(f"/lookup/id/{t['id']}?expand=1")
        filas.append({
            "id": t["id"],
            "nombre": t.get("external_name", ""),
            "exones": len(detalle.get("Exon", [])),
            "aa": (detalle.get("Translation") or {}).get("length", 0),
            "canonico": bool(t.get("is_canonical")),
            "mane": "MANE_Select" in (t.get("tag") or []),
        })
        time.sleep(0.1)   # la API admite ~15 peticiones/s
    imprimir_tabla(filas)


def imprimir_tabla(filas):
    print(f"{'Transcrito':<18}{'Nombre':<12}{'Exones':>7}{'aa':>6}  Notas")
    for f in sorted(filas, key=lambda f: (-f["aa"], f["id"])):
        notas = ", ".join(n for n, v in (("canónico", f["canonico"]), ("MANE Select", f["mane"])) if v)
        print(f"{f['id']:<18}{f['nombre']:<12}{f['exones']:>7}{f['aa']:>6}  {notas}")
    longitudes = sorted({f["aa"] for f in filas if f["aa"]})
    if longitudes:
        print(f"Longitudes de proteína distintas: {len(longitudes)} "
              f"(de {longitudes[0]} a {longitudes[-1]} aa)")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ensembl", metavar="GEN", help="símbolo de un gen humano, p. ej. FGFR2")
    a = p.parse_args()
    if a.ensembl:
        try:
            ensembl(a.ensembl)
        except (urllib.error.HTTPError, urllib.error.URLError) as e:
            raise SystemExit(f"No se pudo consultar Ensembl ({e}). Prueba de nuevo en unos minutos.")
    else:
        gen_juguete()


if __name__ == "__main__":
    main()
