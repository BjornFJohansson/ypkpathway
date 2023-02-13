import gspread
import re

title = "LGM"
name = "MEC_freezer_list"

gc = gspread.service_account()

sheet = gc.open(title)

print(sheet.title)

wks = sheet.worksheet(name)

lol = wks.get_all_values()

regex = r"(?:(?:pYPK(?:a|0))|pTA\d)[^ ]+"

harvest = []
for line in lol:
    for cell in line:
        harvest.extend(re.findall(regex, cell))
        
harvest.sort()


