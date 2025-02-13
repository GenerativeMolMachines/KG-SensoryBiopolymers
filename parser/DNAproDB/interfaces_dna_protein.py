import json
import csv

input_file_path = 'dna-protein.json'  # Путь к вашему входному JSON файлу
output_file_path = 'dna_protein_interfaces.csv'  # Путь к выходному CSV файлу

# Открываем выходной CSV файл для записи
with open(output_file_path, 'w', newline='', encoding='utf-8') as output_file:
    csv_writer = csv.writer(output_file)
    csv_writer.writerow(['structure_id', 'protein_name', 'binding_sites'])  # Записываем заголовки столбцов

    # Открываем входной JSON файл для построчного чтения
    for line in open(input_file_path, 'r', encoding='utf-8'):
        try:
            data = json.loads(line)  # Загружаем каждую строку как JSON объект
            
            # Извлекаем structure_id
            structure_id = data.get("structure_id")  # Получаем structure_id
            
            # Извлекаем protein_metadata
            protein_metadata = data.get("protein_metadata", {})
            for protein_id, metadata in protein_metadata.items():
                protein_name = metadata.get("protein_name")  # Извлекаем protein_name из метаданных
                
                if protein_name is not None:
                    # Извлекаем binding sites из interfaces
                    binding_sites = []  # Список для хранения binding sites
                    interfaces = data.get("interfaces", {}).get("models", [])
                    
                    for model in interfaces:
                        for binding in model:
                            binding_site1 = binding.get("binding_site1")
                            binding_site2 = binding.get("binding_site2")
                            
                            if binding_site1 is not None:
                                binding_sites.append(binding_site1)
                            if binding_site2 is not None:
                                binding_sites.append(binding_site2)

                    # Объединяем все binding sites в одну строку через запятую
                    combined_binding_sites = ', '.join(binding_sites)
                    
                    # Записываем строку в CSV файл
                    csv_writer.writerow([structure_id, protein_name, combined_binding_sites])
                
        except json.JSONDecodeError:
            print("Ошибка при декодировании JSON строки:", line)

print("Данные были успешно сохранены в protein_name.csv.")

