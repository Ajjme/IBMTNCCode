import DownloadDataFromWeatherCompany
import os
import requests
import json
import pandas as pd
import numpy as np
import sys
import decimal
from datetime import datetime as dt, timedelta
import argparse
from pathlib import Path
import time

def getTWCfct_new(latitudes, longitudes, date_begin, date_end):
    print(date_begin, date_end)
    # Number of forecast hours to download for each point
    twc_date_begin = date_begin.strftime("%Y%m%d%H%M")
    twc_date_end = date_end.strftime("%Y%m%d%H%M")
    api_key = 'be01c34853ca4cae81c34853cabcaef1'
    # we need to convert the date_begin/end to correct string
    pointType = 'weighted'
    geocode = '{},{}'.format(latitudes, longitudes)

    url_hod = 'https://api.weather.com/v3/wx/hod/reanalysis/historical/point?pointType={}&geocode={}' \
              '&distance=50&startDateTime={}&endDateTime={}&units=m&format=json&apiKey={}'.format(
        pointType, geocode, twc_date_begin, twc_date_end, api_key)
    print(url_hod)
    r = requests.get(url_hod)
    try:
        weather = json.loads(r.text)
        json.dumps(weather, indent=1)
        df = pd.DataFrame.from_dict(weather)
        df.set_index('observationTimeUtcIso', inplace=True)
        df = df.drop_duplicates()  # some dates are duplicated
        return df
    except ValueError:  # includes simplejson.decoder.JSONDecodeError
        print('Json decoding has failed for', twc_date_begin, twc_date_end)
        return -999, -999

def main_different_observables(start_date, end_date, longitude, latitude):
    date_begin = start_date # since this change inside loop
    year_string = date_begin.strftime("%Y")
    # Create WeatherData folder if doesn't exist
    Path("./WeatherData").mkdir(parents=True, exist_ok=True)
    fname_out = 'WeatherData/TWC_historical_pandas_Full_lat' + str(latitude) \
              + '_lon' + str(longitude)+ "_year" + year_string + '.csv'
    init = True # create or append to pandas data frame
    # We simply iterate over the given date period in 30 day intervals and request data
    # from API. We append to dataframe and end we write to file
    while date_begin < end_date:
        date_end = date_begin + timedelta(days=30)
        if date_end > end_date:   # check that we don't beyond user specified
            date_end = end_date
        # request data from API for given lat, lon, and particular date period
        df_return = getTWCfct_new(latitude, longitude, date_begin, date_end)
        if init:
            df_weather = df_return
            init = False
        else:
            df_weather = df_weather.append(df_return, ignore_index=False)
        date_begin = date_begin + timedelta(days=30)  # TWC allows to download historical
    df_weather.to_csv(fname_out)

        #def float_range(start, stop, step):
        #while start < stop:
           # yield float(start)
           # start += decimal.Decimal(step)

def main():
    #os.system("rm -rf WeatherData/")
    sd = "20150201"
    ed = "20200201"
    start_date = dt.strptime(sd,'%Y%m%d')
    end_date = dt.strptime(ed,'%Y%m%d')

    x = [0,-.1667,-.333, -.667, -.833]
    y = [36.667, 36.83, 37.16, 37.33]
  #x = list(float_range(-2, 0, '0.5'))
  #y = list(float_range(36, 39, '0.5'))

    for num1 in x:
        for num2 in y:
            start_time = time.time()
            print("lat:\t", num1)
            print("long:\t", num2)
            main_different_observables(start_date, end_date, num2, num1)
            print("elapsed time:\n", time.time() - start_time)

if __name__ == '__main__':
    main()



