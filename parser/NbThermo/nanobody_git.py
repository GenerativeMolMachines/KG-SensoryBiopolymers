from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
import time
import os


driver = webdriver.Chrome()


driver.get("https://valdes-tresanco-ms.github.io/NbThermo/")

try:

    wait = WebDriverWait(driver, 30)

    with open("i1.txt", "w", encoding="utf-8") as file:
        table_count = 0

        while table_count < 450:
            clickable_elements = wait.until(
                EC.visibility_of_all_elements_located(
                    (By.CSS_SELECTOR, ".mat-expansion-panel-header")
                )
            )

            if not clickable_elements:

                break

            for element in clickable_elements:
                element.click()

                try:
                    tables = wait.until(
                        EC.presence_of_all_elements_located(
                            (By.CSS_SELECTOR, ".mat-elevation-z1")
                        )
                    )

                    antigen_table = None
                    for table in tables:
                        header = table.find_element(By.TAG_NAME, "h3").text
                        if header == "Antigens":
                            antigen_table = table
                            break

                    if antigen_table:
                        print(element.text)

                        file.write(f"{element.text}" "\n")

                        rows = antigen_table.find_elements(By.TAG_NAME, "tr")
                        all_data = []

                        for row in rows:
                            cells = row.find_elements(By.TAG_NAME, "td")
                            row_data = [
                                cell.text
                                for cell in cells
                                if cell.text.strip()
                                and cell.text.strip() != "open_in_new"
                            ]

                            if row_data:
                                all_data.extend(row_data)

                        line = ", ".join(all_data)
                        print(line)
                        file.write(line + "\n")

                        table_count += 1

                except TimeoutException:
                    print("Таблицы не найдены.")

                time.sleep(1)

finally:

    driver.quit()
