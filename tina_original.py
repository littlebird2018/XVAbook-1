import numpy as np 
from scipy.stats import f 
import pandas as pd 
################################# 
from scipy import stats 
from scipy.stats import f_oneway
import time as time 

df = pd.read_csv(r"D:\Users\data3.csv")
 
def xlMathStat_FstatIncremental(xlValues, xlValueToCheck): 
    try: # Convert input to numpy arrays 
        source_table = np.array(xlValues) 
        value_to_check_range = np.atleast_2d(xlValueToCheck) # Ensure 2D 
        value_to_check = float(value_to_check_range[0, 0]) 
        if np.isnan(value_to_check): 
            raise ValueError("Error: empty field to check") 
        n_row, n_col = source_table.shape 
        if (n_col ==2):
            print ('size issue')
        size_1 = 0 
        size_2 = 0 
        m1_1 = 0.0 
        xi2_1 = 0.0 
        m1_2 = 0.0 
        xi2_2 = 0.0 
        tol = 1.0e-9
        
        for i in range(n_row): 
            for j in range(n_col): 
                field = source_table[i, j] 
                if not np.isnan(field): 
                    vij = float(field) 
                    m1_1 += vij 
                    xi2_1 += vij * vij 
                    if np.abs(vij - value_to_check) > tol: 
                        m1_2 += vij 
                        xi2_2 += vij * vij 
                        size_2 += 1 
                    size_1 += 1 

        if size_2 == size_1: 
            raise ValueError("Value being checked is not in the input range") 
        if size_2 < 2: 
            raise ValueError("Insufficient range size") 
        m1_1 /= size_1 
        m1_2 /= size_2 
        xi2_1 = xi2_1 / size_1 - m1_1 * m1_1 
        xi2_2 = xi2_2 / size_2 - m1_2 * m1_2 
        # Calculate F-statistic 
        df1 = size_1 - 1 
        df2 = size_2 - 1 
        f_stat = (xi2_2 / df2) / (xi2_1 / df1) 
# Calculate p-value 
        p_value = f.cdf(f_stat, df1, df2) 
        return f_stat, p_value 
    except Exception as e: 
        return str(e) 

# Loop through unique values in the 'category' column 
results = {'Tile_Num': [], 'Tile_SequenceNum':[],'Submission': [], 'F_Stat': [], 'P_Value': []}
start=time.time() 
# Loop through unique values in the 'category' column 
results = {'Tile_Num': [], 'Tile_SequenceNum':[],'Submission': [], 'F_Stat': [], 'P_Value': []}
for category in df['Tile_Num'].unique(): # Extract the subset of data for the category 
    category_data = df[df['Tile_Num'] == category] 
    category_values = category_data['Submission'].tolist() 
    category_sequence_num = category_data['Tile_SequenceNum'].tolist() 
    # Loop through the values in the subset 
    for value, sequence_num in zip(category_values, category_sequence_num): 
       
 # Calculate the f-test and p-test for the samples 
        result = xlMathStat_FstatIncremental([category_values], [value]) #append the results to the dictionary 
        results['Tile_Num'].append(category) 
        results['Tile_SequenceNum'].append(sequence_num) # Add Tile_SequenceNum 
        results['Submission'].append(value) 
        results['F_Stat'].append(result[0]) 
        results['P_Value'].append(result[1]) 

results_df = pd.DataFrame(results) 
print (results_df)

end=time.time() 
print("Elapsed (with compilation) = %s"%(end-start)) 
