"""Ejercicio 3: genes del genoma de referencia del SARS-CoV-2 (NC_045512).

Uso:  python ej3_covid19.py

1. Descarga el registro en formato GenBank.
2. Recorre su tabla de características (features).
3. Se queda solo con las de tipo "gene".
4. Lee el nombre del gen en qualifiers.get("gene").
5. Muestra un listado numerado con el nombre y las coordenadas de inicio y fin.

Nota sobre coordenadas: Biopython cuenta desde 0 y el final no se incluye (como en
Python). Para mostrarlas como en GenBank (desde 1 y con el final incluido) se suma 1
al inicio.
"""
from comun import descargar, leer_genbank

ACCESO = "NC_045512"


def genes_del_registro(registro):
    """Lista de (nombre, inicio, fin) de las características de tipo 'gene', en coordenadas GenBank."""
    genes = []
    for caracteristica in registro.features:
        if caracteristica.type == "gene":
            nombre = caracteristica.qualifiers.get("gene", ["(sin nombre)"])[0]
            inicio = int(caracteristica.location.start) + 1
            fin = int(caracteristica.location.end)
            genes.append((nombre, inicio, fin))
    return genes


def main():
    ruta = descargar("nucleotide", ACCESO, "gb", f"{ACCESO}.gb")
    registro = leer_genbank(ruta)
    print(f"\n{registro.id}: {registro.description}")
    print(f"Longitud del genoma: {len(registro.seq)} nucleótidos")
    print(f"Características totales en la tabla: {len(registro.features)}\n")

    genes = genes_del_registro(registro)
    print(f"{'Nº':>3}  {'Gen':<10}{'Inicio':>8}{'Fin':>8}{'Longitud':>10}")
    for i, (nombre, inicio, fin) in enumerate(genes, 1):
        print(f"{i:>3}  {nombre:<10}{inicio:>8}{fin:>8}{fin - inicio + 1:>10}")
    print(f"\nTotal: {len(genes)} genes")


if __name__ == "__main__":
    main()
