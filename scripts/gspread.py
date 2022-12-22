#!/home/bjorn/anaconda3/envs/bjorn38/bin/python
# -*- coding: utf-8 -*-
# ⏳ 😊 😡 🏆 ✅ ❌ ⏱ 👽 😴

from tqdm import tqdm
import gspread
from pydna.utils import cseguid
import datetime
from dateutil import parser
from pandas import DataFrame as df
import pandas

title = "LGM"
name = "pYPKa_A"

gc = gspread.service_account()
sheet = gc.open(title)
wks = sheet.worksheet(name)
lol = wks.get_all_values()
dfrm = df(lol)
