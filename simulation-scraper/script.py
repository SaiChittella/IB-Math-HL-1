from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time

driver = webdriver.Chrome()

try:
    driver.get("https://aa.usno.navy.mil/data/AltAz")


    day = 1

    for day in range(1, 31):
        
        wait = WebDriverWait(driver, 10)

        date_input = wait.until(EC.presence_of_element_located((By.ID, "date")))
        date_input.clear()
        date_input.send_keys(f"07-{day:02d}-2025")

        wait = WebDriverWait(driver, 10)

        # For the interval 
        interval_input = wait.until(EC.presence_of_element_located((By.CLASS_NAME, "intv")))
        interval_input.clear()
        interval_input.send_keys("60")


        # For the latitude and longitude
        latitude_input =  wait.until(EC.presence_of_element_located((By.ID, "lat")))
        latitude_input.clear()
        latitude_input.send_keys("42.60")

        longitude_input = wait.until(EC.presence_of_element_located((By.ID, "lon")))
        longitude_input.clear()
        longitude_input.send_keys("122")

        # For the location
        location_input = wait.until(EC.presence_of_element_located((By.ID, "label")))
        location_input.clear()
        location_input.send_keys('Sammamish')

        # For the timezone
        timezone_input = wait.until(EC.presence_of_element_located((By.ID, "tz")))
        timezone_input.clear()
        timezone_input.send_keys('8')

        # Click on East button
        label = wait.until(EC.element_to_be_clickable((By.XPATH, "//label[text()='East of Greenwich (UTC)']")))
        label.click()

        # Click the "Get Data" button
        get_data_button = driver.find_element(By.ID, "submit")
        get_data_button.click()

        # Wait for the results table to appear (wait for table element)
        tables = wait.until(EC.presence_of_all_elements_located((By.TAG_NAME, "table")))
        table = tables[-1]

        rows = table.find_elements(By.TAG_NAME, "tr")

        for row in rows:
            cells = row.find_elements(By.XPATH, ".//th|.//td")
            row_data = [cell.text.strip() for cell in cells]
            print(row_data)
        
        driver.back()

        day +=1
        print('\n')

finally:
    time.sleep(10) 
    driver.quit()