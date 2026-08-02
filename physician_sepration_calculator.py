#Calculates physician separation based on these criteria:
## When a physician’s NPI was linked to a specific group practice PAC ID (gpid) in year t, but no longer linked to the PAC ID in the following year (t+1), while both the physician (npi) and the group practice (gpid) remained active in the dataset.

import pandas as pd
import os

# 1. Configuration
INPUT_DIR   = 'contains_csvs'
OUTPUT_PATH = ''

# Years represented in the data
ALL_YEARS = list(range(2014, 2024))       # 2014 through 2023
SEP_YEARS = list(range(2014, 2023))       # Only calculate separated for 2014–2022

# Specialty → numeric code
spec_map = {
    'allergy_immunology':        1,
    'cardiology':                2,
    'endocrinology':             3,
    'gastroenterology':          4,
    'geriatric_medicine':        5,
    'hematology_oncology':       6,
    'infectious_disease':        7,
    'general_internal_medicine': 8,
    'nephrology':                9,
    'pulmonary_critical_care':  10,
    'rheumatology':             11
}

# 2. Load and combine specialty CSVs
dfs = []
for fname in os.listdir(INPUT_DIR):
    if fname.endswith('.csv'):
        df = pd.read_csv(os.path.join(INPUT_DIR, fname), dtype=str)
        dfs.append(df)
master = pd.concat(dfs, ignore_index=True)

# 3. Convert data types
master['year'] = master['year'].astype(int)
master['npi']  = master['npi'].astype(str)
master['gpid'] = master['gpid'].astype(str)

# 4. Map specialties to dummy codes
master['specialty'] = master['specialty'].map(spec_map)

# 5. Create year-indexed lookup sets for separated logic
npi_by_year  = { y:set(master.loc[master['year']==y, 'npi'])  for y in ALL_YEARS }
gpid_by_year = { y:set(master.loc[master['year']==y, 'gpid']) for y in ALL_YEARS }
pair_by_year = {
    y:set(zip(
        master.loc[master['year']==y, 'npi'],
        master.loc[master['year']==y, 'gpid']
    ))
    for y in ALL_YEARS
}

# 6. Separation logic: only for years 2014–2022
def flag_sep(row):
    y = row['year']
    if y not in SEP_YEARS:
        return 0  # Do not compute separated beyond 2022
    n, g = row['npi'], row['gpid']
    y_next = y + 1
    return int(
        n in npi_by_year[y_next]
        and g in gpid_by_year[y_next]
        and (n, g) not in pair_by_year[y_next]
    )


master['separated'] = master.apply(flag_sep, axis=1)

# 7. Final cleanup before save
# Remove rows where members is 0 or 1
mask = master['members'].isin([0.0, 1.0, '0', '1', '0.0', '1.0'])
master = master[~mask & master['members'].notna() & (master['members'] != '')]

# Remove rows where gpid is missing or blank
master = master[master['gpid'].notna() & (master['gpid'].str.strip() != '')]

# Remove duplicates
master = master.drop_duplicates(subset=['npi', 'gpid', 'year', 'specialty'])

#Rename npi to NPI
master = master.rename(columns={'npi': 'NPI'})

# 8. Save output
master.to_csv(OUTPUT_PATH, index=False)
print(f"Wrote final dataset with {len(master)} rows to: {OUTPUT_PATH}")
