from bs4 import * 
import requests
import typing as tp
import numpy as np

def make_csv(arr: tp.List[tp.Any], file_path: str) -> bool:
    #try:
        np.savetxt(file_path, np.array(arr).reshape((len(arr[0]), len(arr))), delimiter=',')
        return True
    #except:
        return False

req = requests.get("https://webs.iiitd.edu.in/raghava/prrdb2/browse_sub1.php?token=TLR&col=9")
soup = BeautifulSoup(req.text, "html.parser")
table = soup.find('table')
table_body = table.find('tbody')

data = []
rows = table_body.find_all('tr')
for row in rows:
    cols = row.find_all('td')
    cols = [ele.text.strip() for ele in cols]
    data.append(cols)

f_line = data[0]
open("data.txt", 'w',encoding='utf-8').writelines(f_line)