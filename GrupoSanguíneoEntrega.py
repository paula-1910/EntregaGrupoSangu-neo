from datetime import datetime
import json
import os
import shutil
import matplotlib.pyplot as plt


class GruposangDescriptor:
    """Descriptor para validar que los grupos sanguíneos de los padres sean válidos (A, B, AB, O)."""

    GRUPOS_VALIDOS = {"A", "B", "AB", "O"}

    def __init__(self, nombre_atributo):
        self.nombre_atributo = f"_{nombre_atributo}"

    def __get__(self, instance, owner):
        if instance is None:
            return self
        return getattr(instance, self.nombre_atributo, None)

    def __set__(self, instance, value):
        if not isinstance(value, str):
            raise TypeError("El grupo sanguíneo debe ser una cadena de texto.")

        valor_normalizado = value.strip().upper()
        if valor_normalizado not in self.GRUPOS_VALIDOS:
            raise ValueError(
                f"Grupo sanguíneo '{value}' no válido. Opciones permitidas: {', '.join(sorted(self.GRUPOS_VALIDOS))}"
            )

        setattr(instance, self.nombre_atributo, valor_normalizado)


class CalculadoraHerenciaSangre:
    """Clase principal para la estimación genotípica/fenotípica de grupos sanguíneos,

    persistencia de análisis en JSON y gestión del flujo de trabajo por carpetas.
    """

    padre_grupo = GruposangDescriptor("padre_grupo")
    madre_grupo = GruposangDescriptor("madre_grupo")

    GENOTIPOS_POSIBLES = {
        "A": ["AA", "AO"],
        "B": ["BB", "BO"],
        "AB": ["AB"],
        "O": ["OO"],
    }

    def __init__(self, padre_grupo: str, madre_grupo: str):
        self.padre_grupo = padre_grupo
        self.madre_grupo = madre_grupo
        self._porcentajes = {}

    def _obtener_alelos_posibles(self, grupo: str) -> list[tuple[str, str]]:
        genotipos = self.GENOTIPOS_POSIBLES[grupo]
        return [tuple(g) for g in genotipos]

    def _determinar_fenotipo(self, alelo1: str, alelo2: str) -> str:
        alelos = sorted([alelo1, alelo2])
        genotipo = "".join(alelos)

        if genotipo in ("AA", "AO"):
            return "A"
        elif genotipo in ("BB", "BO"):
            return "B"
        elif genotipo == "AB":
            return "AB"
        elif genotipo == "OO":
            return "O"

    def calcular_probabilidades(self) -> dict[str, float]:
        genotipos_padre = self._obtener_alelos_posibles(self.padre_grupo)
        genotipos_madre = self._obtener_alelos_posibles(self.madre_grupo)

        total_combinaciones = 0
        conteo_fenotipos = {"A": 0, "B": 0, "AB": 0, "O": 0}

        for g_padre in genotipos_padre:
            for g_madre in genotipos_madre:
                for alelo_p in g_padre:
                    for alelo_m in g_madre:
                        fenotipo = self._determinar_fenotipo(alelo_p, alelo_m)
                        conteo_fenotipos[fenotipo] += 1
                        total_combinaciones += 1

        self._porcentajes = {
            grupo: round((conteo / total_combinaciones) * 100, 2)
            for grupo, conteo in conteo_fenotipos.items()
            if conteo > 0
        }

        return self._porcentajes

    def mostrar_resultados(self):
        if not self._porcentajes:
            self.calcular_probabilidades()

        print(
            f"\nResultados para el cruce: Padre ({self.padre_grupo}) x Madre ({self.madre_grupo})"
        )
        print("-" * 55)
        for grupo, porcentaje in self._porcentajes.items():
            print(f"Grupo Sanguíneo {grupo}: {porcentaje}%")
        print("-" * 55)

    def graficar_resultados(self):
        """Muestra un gráfico circular (pie chart) con los porcentajes usando matplotlib."""
        if not self._porcentajes:
            self.calcular_probabilidades()

        etiquetas = [f"Grupo {g}" for g in self._porcentajes.keys()]
        valores = list(self._porcentajes.values())
        colores = ["#d19b26", "#1172D3", "#7aae7a", "#c21dce"]

        plt.figure(figsize=(6, 6))
        plt.pie(
            valores,
            labels=etiquetas,
            autopct="%1.1f%%",
            startangle=140,
            colors=colores[: len(valores)],
            explode=[0.04] * len(valores),
        )
        plt.title(
            f"Probabilidad de Grupo Sanguíneo en la Descendencia\nPadres: {self.padre_grupo} x {self.madre_grupo}"
        )
        plt.tight_layout()
        plt.show()

    # FASE 3: FLUJO COMPLETO DE PROCESAMIENTO DE ARCHIVOS EN DATA

    @classmethod
    def procesar_flujo_archivos(
        cls, base_dir: str = None, mostrar_graficos: bool = False
    ):
        """Procesa archivos .json en 'data/pending', genera los resultados en 'data/results',

        mueve a 'data/done' y elimina de 'data/pending'.
        """
        if base_dir is None:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            base_dir = os.path.join(script_dir, "data")

        pending_dir = os.path.join(base_dir, "pending")
        done_dir = os.path.join(base_dir, "done")
        results_dir = os.path.join(base_dir, "results")

        for carpeta in [pending_dir, done_dir, results_dir]:
            os.makedirs(carpeta, exist_ok=True)

        archivos_pendientes = [
            f for f in os.listdir(pending_dir) if f.endswith(".json")
        ]

        if not archivos_pendientes:
            print(f" No hay archivos pendientes en '{pending_dir}'.")
            return

        print(
            f" Procesando {len(archivos_pendientes)} archivo(s) en '{pending_dir}'..."
        )

        for nombre_archivo in archivos_pendientes:
            ruta_pending = os.path.join(pending_dir, nombre_archivo)
            ruta_done = os.path.join(done_dir, nombre_archivo)
            ruta_results = os.path.join(
                results_dir, f"resultado_{nombre_archivo}"
            )

            try:
                with open(ruta_pending, "r", encoding="utf-8") as f:
                    datos = json.load(f)

                if isinstance(datos, dict):
                    datos = [datos]

                resultados_analisis = []

                for registro in datos:
                    padre = registro.get("padre")
                    madre = registro.get("madre")

                    calculadora = cls(padre_grupo=padre, madre_grupo=madre)
                    porcentajes = calculadora.calcular_probabilidades()

                    if mostrar_graficos:
                        calculadora.graficar_resultados()

                    resultados_analisis.append(
                        {
                            "fecha_analisis": datetime.now().strftime(
                                "%Y-%m-%d %H:%M:%S"
                            ),
                            "padres": {"padre": padre, "madre": madre},
                            "resultados_descendencia": porcentajes,
                        }
                    )

                with open(ruta_results, "w", encoding="utf-8") as f:
                    json.dump(
                        resultados_analisis, f, indent=4, ensure_ascii=False
                    )

                shutil.copy2(ruta_pending, ruta_done)
                os.remove(ruta_pending)

                print(
                    f"Procesado correctamente: '{nombre_archivo}' -> Copiado a 'done/', resultados en 'results/'"
                )

            except Exception as e:
                print(
                    f"Error al procesar el archivo '{nombre_archivo}': {e}"
                )

# BLOQUE PRINCIPAL DE EJECUCIÓN
if __name__ == "__main__":
    #DEMOSTRACIÓN FASE 1: Cálculo e impresión del gráfico
    print("=== DEMOSTRACIÓN FASE 1 ===")
    ejemplo_fase1 = CalculadoraHerenciaSangre(padre_grupo="A", madre_grupo="B")
    ejemplo_fase1.mostrar_resultados()
    ejemplo_fase1.graficar_resultados()

    #DEMOSTRACIÓN FASE 3: Flujo de procesamiento por carpetas
    print("\n=== DEMOSTRACIÓN FASE 3 ===")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    dir_pending = os.path.join(script_dir, "data", "pending")
    os.makedirs(dir_pending, exist_ok=True)

    archivo_demo = os.path.join(dir_pending, "ejemplo_padres.json")
    if not os.path.exists(archivo_demo) and len(os.listdir(dir_pending)) == 0:
        datos_prueba = [
            {"padre": "A", "madre": "O"},
            {"padre": "AB", "madre": "B"},
            {"padre": "O", "madre": "O"},
        ]
        with open(archivo_demo, "w", encoding="utf-8") as f:
            json.dump(datos_prueba, f, indent=4)
        print(
            f"Archivo de prueba creado en '{archivo_demo}' para verificar la ejecución."
        )

    # Procesar archivos de la Fase 3 (pon 'mostrar_graficos=True' si quieres ver el gráfico de cada pareja)
    CalculadoraHerenciaSangre.procesar_flujo_archivos(mostrar_graficos=False)