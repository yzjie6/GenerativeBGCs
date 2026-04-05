import requests
import sys

url = "https://antismash.secondarymetabolites.org/api/v2.0/upload/"
files = {'data': open('results/gbk/TGT_00242.gbk', 'rb')}
data = {'contact': 'developer@project.local'}
try:
    r = requests.post(url, files=files, data=data)
    print(r.status_code, r.text)
except Exception as e:
    print("Error:", e)
