import pandas as pd
import numpy as np
import os

xl_file = "/Volumes/KeithSSD/CB_V4/otu_data/mixing_data/Particle_Data.xlsx"

assert os.path.exists(xl_file)

unformatted_df = pd.read_excel(xl_file, sheet_name='2016')

# the first six columns are date time depth station collection agency 
sample_pts = unformatted_df.iloc[: , :6].dropna()

# the following 300 rows are the point locations 
index_ranges = {i:list(range(i+1,i+301)) for i in sample_pts.index}

# the header for this data is set down into the first column 
second_header = unformatted_df.iloc[0 , 6:].values

# pull out this data and add the right column names 
particle_data = unformatted_df.iloc[1:, 6:].dropna()
particle_data.columns = second_header

assert set([j for i in index_ranges.values() for j in i]) == set(particle_data.index)

# longitude is wrong
particle_data = particle_data.apply(pd.to_numeric, axis=1)
particle_data['Lon'] = particle_data['Lon']/10


one_sample = sample_pts.loc[0, :]
one_section = particle_data.loc[index_ranges[0], :]

 
 np.percentile(, [25, 75])
map = Basemap(llcrnrlon=3.75,llcrnrlat=39.75,urcrnrlon=4.35,urcrnrlat=40.15, resolution = 'h', epsg=5520)
 0.5*((x1*y2 + x2*y3 + x3*y4 + x4*y1) - (x2*y1 + x3*y2 + x4*y3 + x1*y4)) 







