from Bio import SeqIO, AlignIO
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# --- 1. Combinar archivos FASTA en uno solo ---
archivos_fasta = [f"sequence ({i}).fasta" for i in range(5, 31)]
output_fasta = "secuencias_unicas.fasta"

with open(output_fasta, "w") as archivo_salida:
    for archivo in archivos_fasta:
        for secuencia in SeqIO.parse(archivo, "fasta"):
            SeqIO.write(secuencia, archivo_salida, "fasta")

print(f"✅ Todas las secuencias han sido combinadas en '{output_fasta}'.")

# --- 2. Leer y visualizar alineamiento ---
archivo_alineamiento = "_out.250207034227455MdWkmzDeLfXS1HuBU5RgTlsfnormal.aln"
alineamiento = AlignIO.read(archivo_alineamiento, "clustal")

# Mostrar alineamiento en formato gráfico (opcional)
plt.figure(figsize=(12, 6))
plt.text(0, 0.5, "\n".join(str(registro.seq) for registro in alineamiento),
         fontsize=8, family="monospace", wrap=True)
plt.axis('off')
plt.title("Alineamiento de Secuencias")
plt.show()

# --- 3. Cargar y visualizar la matriz de distancias ---
matriz_distancias = pd.read_csv('matriz_distancias.csv', index_col=0)

# Reemplazar comas por puntos y convertir a flotante
matriz_distancias = matriz_distancias.applymap(lambda x: str(x).replace(',', '.')).astype(float)

# Crear el heatmap de distancias genéticas
plt.figure(figsize=(10, 8))
sns.heatmap(matriz_distancias, cmap="YlGnBu", annot=True, fmt=".2f", cbar_kws={'label': 'Distancia Genética'})
plt.title("Matriz de Distancias Genéticas")
plt.xlabel("Secuencias")
plt.ylabel("Secuencias")
plt.show()

# --- 4. Detección de SNPs ---
alignment_matrix = pd.DataFrame([list(str(record.seq)) for record in alineamiento], 
                                index=[record.id for record in alineamiento])

# Identificar posiciones con variaciones (SNPs)
snp_positions = alignment_matrix.nunique(axis=0) > 1
snp_sites = alignment_matrix.loc[:, snp_positions]

# Guardar los SNPs detectados en un CSV
snp_sites.T.to_csv("SNPs_detectados.csv")
print(f"✅ Se encontraron {snp_sites.shape[1]} SNPs. Guardado en 'SNPs_detectados.csv'.")

# --- 5. Visualización de SNPs ---
df = pd.read_csv("SNPs_detectados.csv")

# Convertir valores '-' en NaN para análisis
df.replace('-', pd.NA, inplace=True)

# Contar SNPs por posición
snp_counts = df.iloc[:, 1:].notna().sum(axis=1)

# Gráfico de barras: Frecuencia de SNPs
plt.figure(figsize=(12, 6))
plt.bar(df.iloc[:, 0], snp_counts, color='royalblue')
plt.xlabel("Posición en el Genoma")
plt.ylabel("Frecuencia de SNPs")
plt.title("Frecuencia de SNPs en el Genoma")
plt.show()

# 🔹 **Gráfico de Dispersión: Presencia de SNPs en cada posición**  
plt.figure(figsize=(12, 6))
for i in range(1, df.shape[1]):  # Iterar sobre las secuencias
    posiciones_snp = df.iloc[:, i].notna()
    plt.scatter(df.iloc[:, 0][posiciones_snp], [i] * posiciones_snp.sum(), alpha=0.6)

plt.xlabel("Posición en el Genoma")
plt.ylabel("Número de Secuencia")
plt.title("Distribución de SNPs en Secuencias")
plt.show()
