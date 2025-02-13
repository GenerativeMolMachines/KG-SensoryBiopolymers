import json
import csv

input_file_path = 'dna-protein.json' 
output_file_path = 'dna_protein_interfaces.csv'  


with open(output_file_path, 'w', newline='', encoding='utf-8') as output_file:
    csv_writer = csv.writer(output_file)
    csv_writer.writerow(['structure_id', 'protein_name', 'binding_sites'])  

   
    for line in open(input_file_path, 'r', encoding='utf-8'):
        try:
            data = json.loads(line) 
            
           
            structure_id = data.get("structure_id")  
            
           
            protein_metadata = data.get("protein_metadata", {})
            for protein_id, metadata in protein_metadata.items():
                protein_name = metadata.get("protein_name")  
                
                if protein_name is not None:
                    
                    binding_sites = [] 
                    interfaces = data.get("interfaces", {}).get("models", [])
                    
                    for model in interfaces:
                        for binding in model:
                            binding_site1 = binding.get("binding_site1")
                            binding_site2 = binding.get("binding_site2")
                            
                            if binding_site1 is not None:
                                binding_sites.append(binding_site1)
                            if binding_site2 is not None:
                                binding_sites.append(binding_site2)

                   
                    combined_binding_sites = ', '.join(binding_sites)
                    
                   
                    csv_writer.writerow([structure_id, protein_name, combined_binding_sites])
                
      



