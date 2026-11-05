"""Ejercicio 4: a la caza de proteínas de resistencia en el plásmido pT181 de S. aureus.

Uso:
    python ej4_bacterias.py                      # analiza pT181 (NZ_CP167936.1)
    python ej4_bacterias.py --buscar             # busca pT181 en el NCBI por nombre y organismo
    python ej4_bacterias.py --acceso NC_001371   # analiza el código de acceso que se indique

1. Localiza y descarga del NCBI el plásmido en formato GenBank.
2. Busca todas las regiones codificantes (características de tipo "CDS").
3. Extrae la secuencia de cada CDS (teniendo en cuenta en qué hebra está) y la traduce
   con la tabla de código genético bacteriana (tabla 11).
4. Compara la traducción con la que viene anotada en el propio registro y guarda las
   proteínas en resultados/proteinas_<acceso>.fasta.

Nota: el enunciado da el código NC_001371, pero en el NCBI ese registro es el plásmido
ColE1 de Escherichia coli, no pT181. Por eso, por defecto, se analiza NZ_CP167936.1, el
registro de pT181 de S. aureus que devuelve la búsqueda en el NCBI (se fija para que el
resultado sea siempre el mismo). Con --buscar se repite la búsqueda y con
--acceso NC_001371 se analiza el registro del enunciado.
"""
import argparse
import os

from Bio import Entrez, SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.SeqRecord import SeqRecord

from comun import CARPETA_RESULTADOS, consultar_ncbi, descargar, leer_genbank

ACCESO_ENUNCIADO = "NC_001371"
ACCESO_PT181 = "NZ_CP167936.1"  # pT181 de S. aureus, encontrado con buscar_pt181()
TABLA_BACTERIANA = 11
PALABRAS_RESISTENCIA = ["resist", "tetracycline", "efflux", "lactamase", "chloramphenicol", "erythromycin", "tet"]


def buscar_pt181():
    """Busca en NCBI Nucleotide el plásmido pT181 de S. aureus y devuelve su código de acceso."""
    for termino in ["Staphylococcus aureus[Organism] AND pT181[Title] AND srcdb_refseq[Properties]",
                    "Staphylococcus aureus[Organism] AND pT181[Title] AND complete sequence[Title]"]:
        print(f"Buscando: {termino}")
        handle = consultar_ncbi(Entrez.esearch, db="nucleotide", term=termino)
        ids = Entrez.read(handle)["IdList"]
        handle.close()
        if ids:
            handle = consultar_ncbi(Entrez.esummary, db="nucleotide", id=ids[0])
            resumen = Entrez.read(handle)[0]
            handle.close()
            print(f"Encontrado: {resumen['AccessionVersion']} - {resumen['Title']}\n")
            return resumen["AccessionVersion"]
    raise SystemExit("No se ha encontrado pT181 en el NCBI. Prueba con --acceso <código>.")


def traducir_cds(caracteristica, secuencia_plasmido):
    """Traduce una CDS. extract() da la vuelta a la secuencia si el gen está en la hebra complementaria."""
    adn = caracteristica.extract(secuencia_plasmido)
    tabla = int(caracteristica.qualifiers.get("transl_table", [TABLA_BACTERIANA])[0])
    try:
        # cds=True comprueba que empieza por un codón de inicio y acaba en uno de parada
        return adn.translate(table=tabla, cds=True)
    except TranslationError:
        return adn.translate(table=tabla, to_stop=True)


def posicion(caracteristica):
    """Coordenadas como en GenBank. Si la CDS cruza el origen del plásmido circular, tiene dos tramos."""
    tramos = [f"{int(t.start) + 1}-{int(t.end)}" for t in caracteristica.location.parts]
    return " + ".join(tramos) + (" (cruza el origen)" if len(tramos) > 1 else "")


def main():
    p = argparse.ArgumentParser(description="CDS y proteínas de un plásmido")
    p.add_argument("--acceso", help=f"código de acceso a analizar (el enunciado indica {ACCESO_ENUNCIADO})")
    p.add_argument("--buscar", action="store_true", help="buscar pT181 en el NCBI en vez de usar el código fijado")
    a = p.parse_args()

    if a.acceso:
        acceso = a.acceso
    elif a.buscar:
        acceso = buscar_pt181()
    else:
        acceso = ACCESO_PT181
    ruta = descargar("nucleotide", acceso, "gb", f"{acceso}.gb")
    plasmido = leer_genbank(ruta)
    organismo = plasmido.annotations.get("organism", "?")
    print(f"\n{plasmido.id}: {plasmido.description}")
    print(f"Organismo: {organismo} | longitud: {len(plasmido.seq)} pb | "
          f"topología: {plasmido.annotations.get('topology', '?')}\n")
    if "Staphylococcus" not in organismo:
        print("AVISO: este registro no es de Staphylococcus aureus, así que no es el plásmido pT181.\n")

    proteinas = []
    cds = [c for c in plasmido.features if c.type == "CDS"]
    print(f"Regiones codificantes (CDS) encontradas: {len(cds)}\n")
    for i, c in enumerate(cds, 1):
        gen = c.qualifiers.get("gene", ["-"])[0]
        producto = c.qualifiers.get("product", ["(sin anotar)"])[0]
        id_proteina = c.qualifiers.get("protein_id", [f"CDS_{i}"])[0]
        hebra = "+" if c.location.strand == 1 else "-"
        proteina = traducir_cds(c, plasmido.seq)
        anotada = c.qualifiers.get("translation", [None])[0]

        print(f"CDS {i}: gen {gen} | {producto}")
        print(f"  Posición {posicion(c)} (hebra {hebra}) | {len(proteina)} aminoácidos | proteína {id_proteina}")
        if anotada is not None:
            print(f"  ¿Coincide con la traducción anotada en GenBank? {str(proteina) == anotada}")
        print(f"  {proteina[:60]}{'...' if len(proteina) > 60 else ''}\n")
        proteinas.append(SeqRecord(proteina, id=id_proteina, description=f"{gen} {producto} [{acceso}]"))

    salida = os.path.join(CARPETA_RESULTADOS, f"proteinas_{acceso}.fasta")
    SeqIO.write(proteinas, salida, "fasta")
    print(f"Proteínas guardadas en {salida}")

    resistencia = [r for r in proteinas if any(k in r.description.lower() for k in PALABRAS_RESISTENCIA)]
    print("\nPosibles proteínas de resistencia a antibióticos:")
    for r in resistencia:
        print(f"  - {r.id}: {r.description}")
    if not resistencia:
        print("  Ninguna: ninguna CDS está anotada como proteína de resistencia.")


if __name__ == "__main__":
    main()
