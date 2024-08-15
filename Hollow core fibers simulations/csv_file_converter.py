# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 14:05:26 2024

@author: Owner
"""

import csv

# Read the content of the text file
with open('Glosseries and Meals in French.txt', 'r', encoding='utf-8') as file:
    lines = file.readlines()

# Prepare data for CSV
data = []
for line in lines:
    # Extract the French word and its English translation
    french_word = line.split('(')[0].strip()
    english_translation = line.split('(')[1].split(')')[0].split(',')[0].strip()
    data.append([french_word, english_translation])

# Write data to a CSV file
with open('output.csv', 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Front', 'Back'])
    writer.writerows(data)

print("CSV file has been created successfully.")