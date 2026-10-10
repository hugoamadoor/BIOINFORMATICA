"""Traducción (ARN -> proteína).

Modelo bacteriano simplificado:
  - La subunidad pequeña del RIBOSOMA se une al ARNm en la secuencia Shine-Dalgarno
    (AGGAGG) y busca el primer codón de inicio AUG que hay después.
  - El ARNt iniciador, cargado con formil-metionina (fMet), se coloca en el sitio P.
  - En cada ciclo de elongación:
      1. entra en el sitio A el ARNt cuyo ANTICODÓN es complementario al codón;
      2. el ribosoma (peptidil transferasa) une la cadena al nuevo aminoácido;
      3. el ribosoma avanza un codón (translocación): el ARNt vacío sale por el sitio E.
  - Al llegar a un codón de paro (UAA, UAG, UGA) no entra ningún ARNt: entra un factor
    de liberación, se suelta la proteína y el ribosoma se separa.
"""
from . import pantalla
from .codigo_genetico import CODIGO_GENETICO, FACTOR_LIBERACION, UNA_LETRA

SHINE_DALGARNO = "AGGAGG"
COMPLEMENTO_ARN_ARN = {"A": "U", "U": "A", "G": "C", "C": "G"}


def anticodon(codon):
    """Anticodón del ARNt, escrito 3'->5' para que quede enfrente del codón."""
    return "".join(COMPLEMENTO_ARN_ARN[base] for base in codon)


def buscar_inicio(arnm):
    """Devuelve (posición de Shine-Dalgarno o -1, posición del primer AUG después)."""
    sd = arnm.find(SHINE_DALGARNO)
    desde = sd + len(SHINE_DALGARNO) if sd != -1 else 0
    return sd, arnm.find("AUG", desde)


def dibujar_ribosoma(codones, a):
    """Dibuja los codones alrededor del ribosoma, marcando los sitios E, P y A."""
    desde = max(0, a - 3)
    hasta = min(len(codones), a + 4)
    linea = "   ARNm 5'-" + ("... " if desde > 0 else "")
    sitios = " " * len(linea)
    for i in range(desde, hasta):
        nombre = {a - 2: "E", a - 1: "P", a: "A"}.get(i, " ")
        texto = f"[{codones[i]}]" if nombre != " " else f" {codones[i]} "
        linea += pantalla.color(texto, "cian") if nombre != " " else texto
        sitios += f"  {nombre}  "
    linea += ("..." if hasta < len(codones) else "") + "-3'"
    print(linea)
    print(pantalla.color(sitios, "cian"))


def simular_traduccion(arnm):
    """Simula la traducción del ARNm y devuelve la proteína (lista de aminoácidos)."""
    pantalla.titulo("3. TRADUCCIÓN  (ARN -> PROTEÍNA)")
    sd, inicio = buscar_inicio(arnm)
    if inicio == -1:
        print("   No hay ningún codón AUG: el ribosoma no puede empezar la traducción.")
        pantalla.pausa()
        return {"proteina": [], "inicio": -1, "paro": None, "sd": sd}

    if sd != -1:
        print(f"   La subunidad pequeña del ribosoma se une a la secuencia Shine-Dalgarno "
              f"(AGGAGG, posición {sd + 1} del ARNm)\n   y busca el primer AUG que hay después: "
              f"posición {inicio + 1} del ARNm.")
    else:
        print("   No hay secuencia Shine-Dalgarno: el ribosoma empieza en el primer AUG del ARNm,"
              f" en la posición {inicio + 1}.")

    codones = [arnm[i:i + 3] for i in range(inicio, len(arnm) - 2, 3)]
    proteina = []
    paro = None

    for numero, codon in enumerate(codones):
        aminoacido = CODIGO_GENETICO[codon]

        if aminoacido == "STOP":
            paro = codon
            if pantalla.mostrar_detalle():
                pantalla.subtitulo(f"Codón {numero + 1}: {codon} es un codón de PARO")
                dibujar_ribosoma(codones, numero)
                print(f"   No hay ARNt para {codon}: entra el factor de liberación ({FACTOR_LIBERACION[codon]}).")
                print("   Se suelta la proteína y las dos subunidades del ribosoma se separan.")
                pantalla.pausa()
            break

        proteina.append(aminoacido)
        if pantalla.mostrar_detalle():
            if numero == 0:
                pantalla.subtitulo(f"Iniciación: codón {codon}")
                dibujar_ribosoma(codones, 1)   # el AUG queda en el sitio P
                print(f"   El ARNt iniciador (anticodón 3'-{anticodon(codon)}-5') trae formil-metionina"
                      " al sitio P\n   y se une la subunidad grande del ribosoma. El sitio A queda libre,"
                      f" frente al codón {codones[1] if len(codones) > 1 else '-'}.")
            else:
                pantalla.subtitulo(f"Elongación: codón {numero + 1} = {codon} -> {aminoacido}")
                dibujar_ribosoma(codones, numero)
                print(f"   1. Entra en el sitio A el ARNt-{aminoacido} con anticodón 3'-{anticodon(codon)}-5'"
                      f" (empareja con 5'-{codon}-3').")
                print("   2. La peptidil transferasa une la cadena al nuevo aminoácido.")
                print("   3. El ribosoma avanza un codón; el ARNt vacío sale por el sitio E.")
            print("   Cadena: " + pantalla.color("-".join(proteina), "verde"))
            pantalla.pausa()

    if paro is None:
        print("   Atención: se acabó el ARNm sin encontrar un codón de paro.")

    pantalla.subtitulo("Resultado")
    tres_letras = "-".join(proteina) if len(proteina) <= 20 else "-".join(proteina[:20]) + "-..."
    print(f"   Proteína ({len(proteina)} aminoácidos): {tres_letras}")
    print(f"   En código de una letra: {''.join(UNA_LETRA[aa] for aa in proteina)}")
    print("   (en bacterias, la formil-metionina inicial suele perder el grupo formilo después)")
    pantalla.pausa()

    return {"proteina": proteina, "inicio": inicio, "paro": paro, "sd": sd}
