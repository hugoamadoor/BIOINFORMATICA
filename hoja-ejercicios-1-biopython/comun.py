"""Funciones compartidas por los ejercicios: conexión con el NCBI y descarga de registros."""
import os
import time
from urllib.error import HTTPError, URLError

from Bio import Entrez, SeqIO


Entrez.email = os.environ.get("NCBI_EMAIL", "hugo.amador101@alu.ulpgc.es")

CARPETA_RESULTADOS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultados")
os.makedirs(CARPETA_RESULTADOS, exist_ok=True)


def usar_certificados_del_sistema():
    """Evita errores de certificado HTTPS en redes con antivirus o proxy (si truststore está instalado)."""
    try:
        import truststore
        truststore.inject_into_ssl()
    except ImportError:
        pass


usar_certificados_del_sistema()


def consultar_ncbi(funcion, intentos=4, espera=5, **parametros):
    """Llama a una función de Entrez (esearch, esummary, efetch) y reintenta si el NCBI falla.

    El NCBI devuelve a veces errores temporales (500, 502, 503) aunque la consulta sea
    correcta. Biopython reintenta, pero sin esperar entre intentos, así que aquí se espera
    unos segundos (5, 10, 15...) antes de volver a probar.
    """
    for intento in range(1, intentos + 1):
        try:
            return funcion(**parametros)
        except (HTTPError, URLError) as error:
            # Un error 4xx (salvo 429, "demasiadas peticiones") indica una consulta mal hecha: no se reintenta
            if isinstance(error, HTTPError) and 400 <= error.code < 500 and error.code != 429:
                raise
            if intento == intentos:
                raise SystemExit(f"El NCBI no responde ({error}). Es un fallo temporal de su servidor: "
                                 "espera un par de minutos y vuelve a ejecutar el script.")
            print(f"  El NCBI ha devuelto un error ({error}). Reintento {intento} de {intentos - 1} "
                  f"en {espera * intento} s...")
            time.sleep(espera * intento)


def descargar(db, identificador, formato, nombre_fichero):
    """Descarga un registro del NCBI con Entrez.efetch y lo guarda en resultados/.

    Si el fichero ya existe, no se vuelve a descargar.
    db: "nucleotide" o "protein"; formato: "fasta" o "gb".
    """
    ruta = os.path.join(CARPETA_RESULTADOS, nombre_fichero)
    if os.path.exists(ruta):
        print(f"(Usando la copia guardada en {ruta})")
        return ruta
    print(f"Descargando {identificador} de la base de datos '{db}' en formato {formato}...")
    # "gbwithparts" garantiza que el GenBank incluya la secuencia completa: con "gb" algunos
    # registros RefSeq llegan solo con las anotaciones y la secuencia queda indefinida.
    tipo = "gbwithparts" if formato == "gb" else formato
    handle = consultar_ncbi(Entrez.efetch, db=db, id=identificador, rettype=tipo, retmode="text")
    texto = handle.read()
    handle.close()
    with open(ruta, "w", encoding="utf-8") as f:
        f.write(texto)
    return ruta


def leer_genbank(ruta):
    """Lee un GenBank; si se guardó sin secuencia (descarga antigua), lo borra y avisa."""
    registro = SeqIO.read(ruta, "genbank")
    try:
        bytes(registro.seq[:1])
    except Exception:
        os.remove(ruta)
        raise SystemExit(f"{ruta} no contenía la secuencia y se ha borrado. Vuelve a ejecutar el script.")
    return registro
