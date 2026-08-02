import pandas as pd
import os

# --------------------
# 1. Configuration
# --------------------
#Care compare database files (2014-2023) should go into RAW_DIR. All columns in these files should be named as in step 10 below
#HCPSC files (from Physician & Other Practitioners by Provider and Service post-processed by hcpsc_code) should go into the SRVC_DIR. The Rndr_NPI column should be renamed to NPI
#PROCESSED_DIR is the output directory

YEARS = list(range(2014, 2024))
RAW_DIR = r""
SRVC_DIR = r""
PROCESSED_DIR = r""
os.makedirs(PROCESSED_DIR, exist_ok=True)

# Required training duration by specialty
training_years = {
    'general_internal_medicine': 3,
    'allergy_immunology': 5,
    'cardiology': 6,
    'endocrinology': 5,
    'gastroenterology': 6,
    'geriatric_medicine': 4,
    'hematology_oncology': 6,
    'infectious_disease': 5,
    'nephrology': 5,
    'pulmonary_critical_care': 6,
    'rheumatology': 5
}

def strip_all_whitespace(df):
    # for every object (string) column, strip leading/trailing spaces
    for col in df.select_dtypes(include=['object']):
        df[col] = df[col].str.strip()
    return df

# --------------------
# 2. Specialty Filters
# --------------------
def filter_subspecialty(df, subspec):
    sec_cols = ['sec_spec_1', 'sec_spec_2', 'sec_spec_3', 'sec_spec_4']
    internal_med_condition = (
        (df['pri_spec'] == 'INTERNAL MEDICINE') &
        df[sec_cols].eq(subspec).any(axis=1)
    )
    subspec_condition = (
        (df['pri_spec'] == subspec) &
        ~(
            df[sec_cols].eq('PEDIATRIC MEDICINE').any(axis=1) &
            ~df[sec_cols].eq('INTERNAL MEDICINE').any(axis=1)
        )
    )
    return df[internal_med_condition | subspec_condition]

def filter_general_im(df):
    sec_cols = ['sec_spec_1', 'sec_spec_2', 'sec_spec_3', 'sec_spec_4']
    df_gim = df[df['pri_spec'] == 'INTERNAL MEDICINE'].copy()
    sec = df_gim[sec_cols].fillna('')
    extra = sec.apply(lambda row: any(val and val != 'HOSPITAL MEDICINE' for val in row), axis=1)
    ped = sec.eq('PEDIATRIC MEDICINE').any(axis=1)
    return df_gim[~extra & ~ped]

def filter_cardiology(df):
    cardio_specs = [
        'ADULT CONGENITAL HEART DISEASE (ACHD)',
        'ADVANCED HEART FAILURE AND TRANSPLANT CARDIOLOGY',
        'CARDIAC ELECTROPHYSIOLOGY',
        'CARDIOVASCULAR DISEASE (CARDIOLOGY)',
        'INTERVENTIONAL CARDIOLOGY'
    ]
    sec_cols = ['sec_spec_1', 'sec_spec_2', 'sec_spec_3', 'sec_spec_4']
    df = df[df['pri_spec'].isin(cardio_specs + ['INTERNAL MEDICINE'])].copy()
    im_mask = (
        (df['pri_spec'] == 'INTERNAL MEDICINE') &
        df[sec_cols].isin(cardio_specs).any(axis=1)
    )
    df = df[(df['pri_spec'] != 'INTERNAL MEDICINE') | im_mask]
    ped = df[sec_cols].eq('PEDIATRIC MEDICINE').any(axis=1)
    no_im = ~df[sec_cols].eq('INTERNAL MEDICINE').any(axis=1)
    return df[~(ped & no_im)]

def filter_heme_onc(df):
    heme_specs = ['HEMATOLOGY/ONCOLOGY', 'MEDICAL ONCOLOGY', 'HEMATOLOGY']
    sec_cols = ['sec_spec_1', 'sec_spec_2', 'sec_spec_3', 'sec_spec_4']
    df = df[df['pri_spec'].isin(heme_specs + ['INTERNAL MEDICINE'])].copy()
    im_mask = (
        (df['pri_spec'] == 'INTERNAL MEDICINE') &
        df[sec_cols].isin(heme_specs).any(axis=1)
    )
    df = df[(df['pri_spec'] != 'INTERNAL MEDICINE') | im_mask]
    ped = df[sec_cols].eq('PEDIATRIC MEDICINE').any(axis=1)
    no_im = ~df[sec_cols].eq('INTERNAL MEDICINE').any(axis=1)
    return df[~(ped & no_im)]

def filter_pulm_crit(df):
    sec_cols = ['sec_spec_1', 'sec_spec_2', 'sec_spec_3', 'sec_spec_4']
    df = df[df['pri_spec'].isin(['PULMONARY DISEASE', 'INTERNAL MEDICINE', 'CRITICAL CARE (INTENSIVISTS)'])].copy()
    crit_mask = (
        (df['pri_spec'] == 'CRITICAL CARE (INTENSIVISTS)') &
        df[sec_cols].isin(['INTERNAL MEDICINE', 'PULMONARY DISEASE']).any(axis=1)
    )
    df = df[(df['pri_spec'] != 'CRITICAL CARE (INTENSIVISTS)') | crit_mask]
    im_mask = (
        (df['pri_spec'] == 'INTERNAL MEDICINE') &
        df[sec_cols].isin(['CRITICAL CARE (INTENSIVISTS)', 'PULMONARY DISEASE']).any(axis=1)
    )
    df = df[(df['pri_spec'] != 'INTERNAL MEDICINE') | im_mask]
    ped = df[sec_cols].eq('PEDIATRIC MEDICINE').any(axis=1)
    no_im = ~df[sec_cols].eq('INTERNAL MEDICINE').any(axis=1)
    return df[~(ped & no_im)]

# --------------------
# 3. Cohort Containers
# --------------------
subspecialties = [
    'ALLERGY/IMMUNOLOGY', 'ENDOCRINOLOGY', 'GASTROENTEROLOGY',
    'GERIATRIC MEDICINE', 'INFECTIOUS DISEASE', 'NEPHROLOGY', 'RHEUMATOLOGY'
]

cohorts = {
    'general_internal_medicine': [],
    'cardiology': [],
    'hematology_oncology': [],
    'pulmonary_critical_care': [],
}
for ss in subspecialties:
    cohorts[ss.lower().replace('/', '_').replace(' ', '_')] = []

# --------------------
# 4. Load and Filter Each Physician Compare File by Year
# --------------------
for year in YEARS:
    # Load physician compare file
    path = os.path.join(RAW_DIR, f"{year}.csv")
    df = pd.read_csv(path, dtype=str,encoding='latin-1',engine='python',on_bad_lines='skip')
    df = strip_all_whitespace(df)
    df.columns = [col.lower() for col in df.columns]

    # Columns needed
    cols = [
        'npi', 'gndr', 'cred', 'grd_yr',
        'pri_spec', 'sec_spec_1', 'sec_spec_2', 'sec_spec_3', 'sec_spec_4',
        'org_nm', 'org_pac_id', 'num_org_mem', 'st', 'zip'
    ]
    df = df[cols]
    df['year'] = year

# --------------------
# 5. Remove non-physician credentials and non-US locations
# --------------------
    excluded_creds = [
        'NP', 'PA', 'PT', 'CNS', 'CNA', 'CSW', 'DC', 'OD', 'CP', 'CNM',
        'DPM', 'AU', 'MNT', 'OT', 'AA', 'DDS', 'PSY', 'DDM', 'SCW'
    ]
    df = df[~df['cred'].fillna('').str.upper().isin(excluded_creds)]
    df = df[~df['st'].isin(['GU', 'PR', 'VI', 'MP'])]


# --------------------
# 6. Load, Merge, and Filter Provider Service Counts for E&M Services <50
# --------------------
    #Read totalservices_year file from Medicare Physician & Other Practitioners by Provider and Service processing code (hcpsc_code)
    df_srv = pd.read_csv(os.path.join(SRVC_DIR, f"totalservices_{year}.csv"), dtype=str)
    df_srv = strip_all_whitespace(df_srv)
    df_srv.columns = [col.lower() for col in df_srv.columns]
    df_srv['total_tot_srvcs'] = pd.to_numeric(df_srv['total_tot_srvcs'], errors='coerce')
    #Merge Physician Compare and Physician & Other Practitioners by Provider and Service by NPI
    df = df.merge(df_srv[['npi', 'total_tot_srvcs']], how='left', on='npi')
    df['total_tot_srvcs'] = df['total_tot_srvcs'].fillna(0)
    #Filter out NPIs with less than 50 services
    df = df[df['total_tot_srvcs'] >= 50]

# --------------------
# 7. Load inpatientcodes_year file from hcpsc_code, merge, and remove inpatient physicians
# --------------------

    # Load and merge inpatientcodes_year files
    df_inp = pd.read_csv(os.path.join(SRVC_DIR, f"inpatientcodes_{year}.csv"), dtype=str)
    df_inp.columns = [col.lower() for col in df_inp.columns]
    df_inp['majority_inpatient'] = pd.to_numeric(df_inp['majority_inpatient'], errors='coerce')
    #Merge inpatientcodes_year file with physician compare file
    df = df.merge(df_inp[['npi', 'majority_inpatient']], how='left', on='npi')
    #Remove inpatient physicians
    df = df[(df['majority_inpatient'].isna()) | (df['majority_inpatient'] != 1)]

# --------------------
# 8. Apply Specialty Filters created in step 2
# --------------------

    cohorts['general_internal_medicine'].append(filter_general_im(df))
    cohorts['cardiology'].append(filter_cardiology(df))
    cohorts['hematology_oncology'].append(filter_heme_onc(df))
    cohorts['pulmonary_critical_care'].append(filter_pulm_crit(df))
    for ss in subspecialties:
        key = ss.lower().replace('/', '_').replace(' ', '_')
        cohorts[key].append(filter_subspecialty(df, ss))

# --------------------
# 9. Apply Training-Year Filter to Remove Those Still in Training
# --------------------

for name, frames in cohorts.items():
    df_out = pd.concat(frames, ignore_index=True)

    # Filter valid graduation years
    df_out = df_out[df_out['grd_yr'].fillna('').str.isdigit()]
    df_out['grd_yr'] = df_out['grd_yr'].astype(int)
    df_out['year'] = df_out['year'].astype(int)

    # Apply training duration filter
    required_years = training_years[name]
    df_out = df_out[df_out['year'] - df_out['grd_yr'] >= required_years]
    
# --------------------
# 10. Remove Duplicates and Export
# --------------------

    # Assign specialty column
    df_out['specialty'] = name

    # Clean ZIP codes to 5 digits
    df_out['zip'] = df_out['zip'].str[:5]

    # Rename num_org_mem to members
    df_out.rename(columns={'num_org_mem': 'members'}, inplace=True)

    # Rename org_pac_id to gpid
    df_out.rename(columns={'org_pac_id': 'gpid'}, inplace=True)


    final_cols = [
        'npi', 'gndr', 'cred', 'grd_yr',
        'pri_spec', 'sec_spec_1', 'sec_spec_2', 'sec_spec_3', 'sec_spec_4',
        'org_nm', 'gpid', 'members', 'st', 'zip',
        'total_tot_srvcs', 'year', 'specialty'
    ]

    # Remove duplicates at physician-practice-specialty-year level
    df_out = df_out.drop_duplicates(subset=['npi', 'gpid', 'year', 'specialty'])

    # Select and reorder final columns
    df_out = df_out[final_cols]

    # Save output
    out_path = os.path.join(PROCESSED_DIR, f"{name}.csv")
    df_out.to_csv(out_path, index=False)
    print(f"Wrote {out_path} — {len(df_out)} rows after deduplication")