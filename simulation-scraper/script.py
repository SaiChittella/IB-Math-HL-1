from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
import math

driver = webdriver.Chrome()
all_rows = []

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
        latitude_input.send_keys("47.60")

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

            if len(row_data) == 3 and row_data[0] != 'Time':
                if float(row_data[1][:-1]) < 0 or float(row_data[2][:-1]) < 0:
                    continue
                all_rows.append(row_data)
                print(row_data)
        
        driver.back()

        day +=1
        print('\n')

finally:
    time.sleep(5) 
    driver.quit()


# panel_azimuth_degrees = 180
# panel_tilt_degrees = [0,4,7,9,11,12,15,17.5,19,20,21,22.5]

# for panel_tilt in panel_tilt_degrees:
#     panel_tilt_radians = math.radians(panel_tilt)

#     cos_theta_sum = 0.0

#     for row in all_rows:
#         # Angle of incidence
#         sun_altitude_str = row[1]
#         sun_altitude = float(sun_altitude_str[:-1]) 

#         sun_azimuth_str = row[2]
#         sun_azimuth = float(sun_azimuth_str[:-1])

#         cos_theta_i = (math.sin(math.radians(sun_altitude)) * math.cos(panel_tilt_radians)) + (math.cos(math.radians(sun_altitude)) * math.sin(panel_tilt_radians) * math.cos(math.radians(sun_azimuth - panel_azimuth_degrees)))

#         cos_theta_sum += cos_theta_i

#     print(f'{panel_tilt}° Incidence Sum: {cos_theta_sum}')

#     energy = cos_theta_sum * 1.6 * 0.2
#     print(f'{panel_tilt} Energy Sum: {energy} \n') 