# -*- coding: utf-8 -*-
"""Aptagen.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1ZBPIVVmQyWdyeprASEVTgXuQSYt8BGZ-
"""

import requests
from bs4 import BeautifulSoup
import pandas as pd
import time

BASE_URL = "https://www.aptagen.com/apta-index/"
headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
}
MAX_PAGES = 6  # Количество страниц для парсинга

# Список для хранения данных
aptamers = []

# Функция для парсинга детальной страницы
def parse_detail_page(url):
    response = requests.get(url, headers=headers)
    soup = BeautifulSoup(response.text, "html.parser")

    # Функция для безопасного извлечения данных
    def get_text(label):
        tag = soup.find("label", string=label)
        return tag.find_next("span").text.strip() if tag else "N/A"

    # Находим ID
    id_container = soup.find("h3")  # Находим <h3>
    if id_container:
        id_span = id_container.find("span", attrs={"itemprop": "sku"})  # Ищем <span itemprop="sku">
        aptamer_id = id_span.text.strip() if id_span else "N/A"
    else:
        aptamer_id = "N/A"

    # Извлекаем остальные данные
    target = get_text("Target:")
    category = get_text("Antigen/Target Category:")
    affinity = get_text("Affinity (Kd):")
    binding_buffer = get_text("Binding Conditions/Buffer:")
    binding_temp = get_text("Binding Temp:")
    specificity = get_text("Specificity:")
    comments = get_text("Comments:")
    length = get_text("Length:")
    molecular_weight = get_text("Molecular Weight:")
    extinction_coefficient = get_text("Extinction Coefficient:")

    # Извлекаем последовательность (с учетом вложенных span)
    sequence_wrap = soup.find("span", id="aptamer-sequence-wrap")
    if sequence_wrap:
        sequence_span = sequence_wrap.find("span")
        sequence = sequence_span.text.strip() if sequence_span else "N/A"
    else:
        sequence = "N/A"

    return {
        "ID": aptamer_id,
        "Target": target,
        "Category": category,
        "Affinity (Kd)": affinity,
        "Binding Conditions/Buffer": binding_buffer,
        "Binding Temp": binding_temp,
        "Specificity": specificity,
        "Comments": comments,
        "Sequence": sequence,
        "Length": length,
        "Molecular Weight": molecular_weight,
        "Extinction Coefficient": extinction_coefficient,
        "URL": url
    }

# Цикл по страницам
for page in range(1, MAX_PAGES + 1):
    print(f"Парсим страницу {page}...")

    # Формируем URL
    page_url = f"{BASE_URL}?pageNumber={page}&targetCategory=All&affinityUnits=&affinityMin=&affinityMax=&aptamerChemistry=DNA&sortBy=Length+%28low-high%29&aptazyme=None&pageSize=100&searchQuery="
    response = requests.get(page_url, headers=headers)
    soup = BeautifulSoup(response.text, "html.parser")

    # Находим все аптамеры на странице
    aptamer_blocks = soup.find_all("div", class_="row result")

    for block in aptamer_blocks:
        name_tag = block.find("a")
        url = name_tag["href"] if name_tag else "N/A"

        if url != "N/A":
            print(f"Парсим {url}...")
            aptamer_data = parse_detail_page(url)
            aptamers.append(aptamer_data)

            # Добавляем задержку, чтобы избежать блокировки
            time.sleep(1)

# Сохраняем в CSV
df = pd.DataFrame(aptamers)
df.to_csv("aptamers_data.csv", index=False, encoding="utf-8")
print("Данные сохранены в aptamers_data.csv")