"""Ejercicio 6. Pipeline del dogma central: replicación -> transcripción -> traducción.

Uso:
  python src/ej6_pipeline.py data/INS_ENST00000381330_cds.fasta \
         --referencia data/P01308_INS_HUMAN.fasta --salida resultados

El FASTA debe contener una hebra de ADN escrita 5'->3' (por defecto, la codificante).
Cada paso se informa por pantalla y en <salida>/pipeline.log, y los productos se
guardan en FASTA dentro de <salida>/.
"""
import argparse
import logging
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from dogma import (buscar_orf, complementaria, contenido_gc, leer_fasta,
                   replicar, traducir, transcribir_desde_codificante,
                   transcribir_desde_molde)

log = logging.getLogger("pipeline")
PAREJAS = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}


def configurar_log(salida):
    os.makedirs(salida, exist_ok=True)
    fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s", "%H:%M:%S")
    log.setLevel(logging.INFO)
    for h in (logging.StreamHandler(), logging.FileHandler(os.path.join(salida, "pipeline.log"), "w", "utf-8")):
        h.setFormatter(fmt)
        log.addHandler(h)


def guardar(salida, nombre, secuencia, ident, descripcion):
    ruta = os.path.join(salida, nombre)
    SeqIO.write(SeqRecord(Seq(secuencia), id=ident, description=descripcion), ruta, "fasta")
    log.info("    guardado %s", ruta)


def corto(s, n=30):
    return s if len(s) <= 2 * n else f"{s[:n]}...{s[-n:]}"


def paso_lectura(fasta, hebra):
    log.info("[0/3] LECTURA Y VALIDACIÓN")
    ident, s = leer_fasta(fasta)
    log.info("  Registro %s: %d nt, GC = %.1f %%", ident, len(s), contenido_gc(s))
    codificante = s if hebra == "codificante" else complementaria(s)
    if hebra == "molde":
        log.info("  La entrada es la hebra molde: se obtiene la codificante por complementariedad")
    avisos = [(len(codificante) % 3 != 0, "la longitud no es múltiplo de 3"),
              (not codificante.startswith("ATG"), "no empieza por ATG"),
              (codificante[-3:] not in ("TAA", "TAG", "TGA"), "no termina en codón de paro")]
    for condicion, texto in avisos:
        if condicion:
            log.warning("  La secuencia %s: puede no ser una CDS completa", texto)
    if not any(c for c, _ in avisos):
        log.info("  Parece una CDS completa: ATG ... %s, %d codones", codificante[-3:], len(codificante) // 3)
    return ident, codificante


def paso_replicacion(ident, codificante, salida):
    log.info("[1/3] REPLICACIÓN (semiconservativa)")
    log.info("  Helicasa: separa las dos hebras parentales rompiendo los puentes de hidrógeno")
    r = replicar(codificante)
    log.info("  Hebra parental superior 5'->3': %s", corto(r["parental_superior"]))
    log.info("  Hebra parental inferior 5'->3': %s", corto(r["parental_inferior"]))
    log.info("  Primasa + ADN polimerasa: síntesis 5'->3' sobre cada molde (conductora y retardada)")
    log.info("  Nueva hebra sobre el molde superior: %s", corto(r["nueva_sobre_superior"]))
    log.info("  Nueva hebra sobre el molde inferior: %s", corto(r["nueva_sobre_inferior"]))
    log.info("  Ligasa: sella los fragmentos de Okazaki de la hebra retardada")
    for nombre, (a, b) in (("Hija 1", r["hija_1"]), ("Hija 2", r["hija_2"])):
        ok = sum((x, y) in PAREJAS for x, y in zip(a, reversed(b)))
        log.info("  %s: %d/%d pares Watson-Crick correctos", nombre, ok, len(a))
    identicas = r["hija_1"] == r["hija_2"] == (codificante, complementaria(codificante))
    log.info("  Las dos moléculas hijas son idénticas a la parental: %s", identicas)
    guardar(salida, "1_replicacion_nueva_sobre_superior.fasta", r["nueva_sobre_superior"], ident + "_nueva1",
            "hebra nueva sintetizada sobre el molde superior (5'->3')")
    guardar(salida, "1_replicacion_nueva_sobre_inferior.fasta", r["nueva_sobre_inferior"], ident + "_nueva2",
            "hebra nueva sintetizada sobre el molde inferior (5'->3')")
    return r


def paso_transcripcion(ident, codificante, salida):
    log.info("[2/3] TRANSCRIPCIÓN")
    molde = complementaria(codificante)
    log.info("  Hebra molde (leída 3'->5' por la ARN polimerasa), escrita 5'->3': %s", corto(molde))
    arnm = transcribir_desde_molde(molde)
    log.info("  ARNm 5'->3' (T -> U): %s", corto(arnm))
    coincide = arnm == transcribir_desde_codificante(codificante)
    log.info("  Comprobación: el ARNm coincide con la hebra codificante cambiando T por U: %s", coincide)
    if not coincide:
        log.error("  La transcripción no es coherente")
    guardar(salida, "2_transcripcion_arnm.fasta", arnm, ident + "_ARNm", "ARNm 5'->3'")
    return arnm


def paso_traduccion(ident, arnm, salida, referencia=None):
    log.info("[3/3] TRADUCCIÓN (código genético estándar)")
    inicio, paro = buscar_orf(arnm)
    if inicio is None:
        log.error("  No hay codón AUG: no se puede iniciar la traducción")
        return ""
    log.info("  Codón de inicio AUG en la posición %d", inicio + 1)
    if paro is None:
        log.warning("  No hay codón de paro en fase: la proteína podría estar incompleta")
    else:
        log.info("  Codón de paro %s en la posición %d", arnm[paro:paro + 3], paro + 1)
    proteina = traducir(arnm[inicio:])
    pa = ProteinAnalysis(proteina)
    log.info("  Proteína: %d aa, %.1f kDa, pI %.2f", len(proteina), pa.molecular_weight() / 1000,
             pa.isoelectric_point())
    log.info("  Secuencia: %s", corto(proteina))
    if referencia:
        registro = next(SeqIO.parse(referencia, "fasta"))
        ref_id, ref = registro.id, str(registro.seq)
        iguales = sum(a == b for a, b in zip(proteina, ref))
        log.info("  Comparación con %s: %d/%d posiciones idénticas%s", ref_id, iguales, len(ref),
                 " (coincidencia exacta)" if proteina == ref else "")
    guardar(salida, "3_traduccion_proteina.fasta", proteina, ident + "_proteina", "traducción de la ORF")
    return proteina


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fasta")
    p.add_argument("--hebra", choices=["codificante", "molde"], default="codificante")
    p.add_argument("--referencia", help="FASTA de proteína con la que comparar el resultado")
    p.add_argument("--salida", default="resultados")
    a = p.parse_args()

    configurar_log(a.salida)
    log.info("Inicio del pipeline: %s", a.fasta)
    ident, codificante = paso_lectura(a.fasta, a.hebra)
    paso_replicacion(ident, codificante, a.salida)
    arnm = paso_transcripcion(ident, codificante, a.salida)
    paso_traduccion(ident, arnm, a.salida, a.referencia)
    log.info("Pipeline terminado. Resultados en %s/", a.salida)


if __name__ == "__main__":
    main()
