import os
import pandas as pd

def process_csv_files(folder_path):
    output_folder = os.path.join(folder_path, "processed_tables")
    os.makedirs(output_folder, exist_ok=True)
    
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".csv"):
            file_path = os.path.join(folder_path, file_name)
            df = pd.read_csv(file_path)
            
            if "Results_column" in df.columns:
                df["Consensus_classification"] = df["Results_column"].apply(
                    lambda x: "Conserved_50" if x in ["Conserved", "Core", "Soft core"] else ""
                )
                
                output_path = os.path.join(output_folder, file_name)
                df.to_csv(output_path, index=False)
                print(f"Processed and saved: {output_path}")
            else:
                print(f"Skipping {file_name}: 'Results_column' not found.")

# Ruta de la carpeta con los archivos CSV
folder_path = "C:/Users/gonza/OneDrive/Escritorio/Doctorado/Tablas_todos/Tablas_todos"
process_csv_files(folder_path)