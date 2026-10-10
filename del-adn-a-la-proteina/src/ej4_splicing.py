"""Ejercicio 4. Splicing alternativo.

Uso:  python src/ej4_splicing.py                 # gen de juguete con 5 exones
      python src/ej4_splicing.py --ensembl FGFR2  # isoformas reales (requiere Internet)
      python src/ej4_splicing.py --ensembl FGFR2 --recargar  # ignora la consulta guardada

Parte 1: los exones 1 y 5 son constitutivos y los exones 2, 3 y 4 pueden incluirse
o saltarse, lo que da 2^3 = 8 isoformas. Para cada una se comprueba si conserva la
pauta de lectura, dónde aparece el primer codón de paro y si el ARNm sería
candidato a degradación por NMD (regla de los 50-55 nt).
Parte 2: consulta la API REST de Ensembl y compara los transcritos de un gen humano.
"""
import argparse
import itertools
import json
import os
import ssl
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


MENSAJE_SSL = (
    "No se pudo verificar el certificado HTTPS de Ensembl ({}).\n"
    "Suele deberse a un antivirus o a la red (universidad, empresa) que inspecciona el\n"
    "tráfico HTTPS con su propio certificado. Solución: instala 'truststore'\n"
    "(pip install truststore) para que Python use los certificados del sistema, o\n"
    "ejecuta el script desde otra red.")


def usar_certificados_del_sistema():
    """Hace que Python confíe en los certificados del sistema operativo (Windows/macOS).

    Con 'truststore' instalado se aceptan también las CA que añaden antivirus o redes
    corporativas; sin él, Python usa solo su propio almacén y puede fallar la verificación.
    """
    try:
        import truststore
        truststore.inject_into_ssl()
    except ImportError:
        pass


def pedir_json(ruta, datos=None, reintentos=4, timeout=30):
    """Petición a la API REST de Ensembl (GET, o POST si se pasan datos) con reintentos
    ante errores temporales: 429 (demasiadas peticiones), 5xx, cortes de red y timeouts."""
    url = f"{API}{ruta}"
    cabeceras = {"Content-Type": "application/json", "Accept": "application/json"}
    cuerpo = json.dumps(datos).encode() if datos is not None else None
    for intento in range(1, reintentos + 1):
        try:
            peticion = urllib.request.Request(url, data=cuerpo, headers=cabeceras)
            with urllib.request.urlopen(peticion, timeout=timeout) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503, 504) or intento == reintentos:
                raise
            motivo = f"HTTP {e.code}"
            espera = float(e.headers.get("Retry-After") or 2 * intento)
        except urllib.error.URLError as e:
            if isinstance(e.reason, ssl.SSLCertVerificationError):
                raise SystemExit(MENSAJE_SSL.format(e.reason)) from None
            if intento == reintentos:
                raise
            motivo, espera = f"error de red: {e.reason}", 2 * intento
        except TimeoutError:
            if intento == reintentos:
                raise
            motivo, espera = f"sin respuesta en {timeout} s", 2 * intento
        print(f"    Ensembl: {motivo}; reintento {intento}/{reintentos - 1} en {espera:.0f} s")
        time.sleep(espera)


def detalles_transcritos(ids, lote=10):
    """Exones y longitud de la proteína de cada transcrito.

    Se piden en lotes con POST /lookup/id (una petición por cada 10 transcritos en vez
    de una por transcrito). Si un lote falla, se piden sus transcritos uno a uno.
    """
    detalles = {}
    for i in range(0, len(ids), lote):
        grupo = ids[i:i + lote]
        print(f"  Transcritos {i + 1}-{i + len(grupo)} de {len(ids)}...")
        try:
            detalles.update(pedir_json("/lookup/id", {"ids": grupo, "expand": 1}) or {})
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
            print("    El lote ha fallado; se piden uno a uno")
            for t in grupo:
                detalles[t] = pedir_json(f"/lookup/id/{t}?expand=1")
        time.sleep(0.2)
    return detalles


def ensembl(simbolo, recargar=False):
    """Lista los transcritos codificantes de un gen humano con su nº de exones y longitud.

    No se usa lookup/symbol?expand=1 sobre el gen completo porque en genes con muchos
    transcritos (FGFR2 tiene 59) la API devuelve a menudo un error 500. En su lugar:
    1) lookup del gen, 2) overlap para listar sus transcritos, 3) detalles por lotes.
    El resultado se guarda en resultados/ensembl_<GEN>.json para no repetir la consulta.
    """
    cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "resultados",
                         f"ensembl_{simbolo.upper()}.json")
    if os.path.exists(cache) and not recargar:
        with open(cache, encoding="utf-8") as f:
            guardado = json.load(f)
        print(f"Usando la consulta guardada en {os.path.normpath(cache)} "
              f"({guardado['fecha']}). Usa --recargar para repetirla.")
        print(guardado["resumen"])
        imprimir_tabla(guardado["filas"])
        return

    usar_certificados_del_sistema()
    print(f"Consultando Ensembl ({API}) para {simbolo}...")
    gen = pedir_json(f"/lookup/symbol/homo_sapiens/{simbolo}")
    transcritos = [t for t in pedir_json(f"/overlap/id/{gen['id']}?feature=transcript")
                   if t.get("Parent") == gen["id"]]
    codif = [t for t in transcritos if t.get("biotype") == "protein_coding"]
    resumen = (f"{gen.get('display_name')} ({gen['id']}), cromosoma {gen.get('seq_region_name')}, "
               f"hebra {'+' if gen.get('strand') == 1 else '-'}: {len(transcritos)} transcritos, "
               f"{len(codif)} codificantes de proteína")
    print(resumen + ". Descargando detalles...")
    detalles = detalles_transcritos([t["id"] for t in codif])
    filas = []
    for t in codif:
        d = detalles.get(t["id"]) or {}
        filas.append({
            "id": t["id"],
            "nombre": t.get("external_name", ""),
            "exones": len(d.get("Exon", [])),
            "aa": (d.get("Translation") or {}).get("length", 0),
            "canonico": bool(t.get("is_canonical")),
            "mane": "MANE_Select" in (t.get("tag") or []),
        })
    imprimir_tabla(filas)
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "w", encoding="utf-8") as f:
        json.dump({"fecha": time.strftime("%Y-%m-%d %H:%M"), "resumen": resumen, "filas": filas},
                  f, ensure_ascii=False, indent=1)
    print(f"Consulta guardada en {os.path.normpath(cache)}")


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
    p.add_argument("--recargar", action="store_true", help="repite la consulta aunque haya una guardada")
    a = p.parse_args()
    if a.ensembl:
        try:
            ensembl(a.ensembl, a.recargar)
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as e:
            raise SystemExit(f"No se pudo consultar Ensembl ({e}). Prueba de nuevo en unos minutos.")
    else:
        gen_juguete()


if __name__ == "__main__":
    main()
