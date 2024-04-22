import pandas as pd
import numpy as np
import pyarrow.parquet as pq
import argparse

def main(args):
    # Extract the number of rows from the original file metadata
    nrows = pq.ParquetFile(args.input_file).metadata.num_rows
    print("Number of rows in original file:", nrows)

    # Read the schema and data of the parquet file
    schema = pq.read_schema(args.input_file, memory_map=True)
    df = pd.read_parquet(args.input_file)

    print("Original df summary ----")
    print(df.head())
    print(df.shape)
    print(df.dtypes)
    print(df.describe())
    
    # Initialize new objects
    numerical_cols = []
    categorical_cols = []
    other_cols = []
    df_fake = pd.DataFrame()

    # Generate fake data
    for field in schema:
        x = df[field.name]
        if field.type in ['float32', 'float64']:
            numerical_cols.append(field.name)
            df_fake[field.name] = np.random.default_rng().uniform(low=min(x), high=max(x), size=args.sample_nrows)
            missing_value = np.nan  # Suitable for numerical data

        elif field.type in ['int32', 'int64']:
            categorical_cols.append(field.name)
            df_fake[field.name] = np.random.choice(x.dropna().unique(), size=args.sample_nrows, replace=True)
            missing_value = pd.NA  # pandas' nullable integer type supports pd.NA

        elif field.type == 'string':
            categorical_cols.append(field.name)
            df_fake[field.name] = np.random.choice(x.dropna().unique(), size=args.sample_nrows, replace=True)
            missing_value = None  # or "" if you prefer to use an empty string to represent missing strings

        else:
            other_cols.append(field.name)
            df_fake[field.name] = np.random.choice(x.dropna().unique(), size=args.sample_nrows, replace=True)
            missing_value = None  # Assuming 'None' as a generic missing value for other types

        # Inject the same proportion of missing values in the original data into the fake data
        if args.inject_missing:
            missing_proportion = x.isnull().sum() / len(x)
            n_missing = int(missing_proportion * args.sample_nrows)
            missing_indices = np.random.choice(df_fake.index, n_missing, replace=False)
            df_fake.loc[missing_indices, field.name] = missing_value

        print(f"Field: {field.name}, Type: {field.type}, Missing Proportion: {missing_proportion}")

    print("Categorical Columns:", categorical_cols)
    print("Numerical Columns:", numerical_cols)
    print("Other Columns:", other_cols)

    print("Fake df summary ----")
    print(df_fake.head())
    print(df_fake.shape)
    print(df_fake.dtypes)
    print(df_fake.describe())

    # Save the fake data to a parquet file
    df_fake.to_parquet(f"{args.output_prefix}.parquet")
    df_fake.to_csv(f"{args.output_prefix}.csv", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", 
                        default = "data/processed/data.parquet", 
                        type= str
                       )
    parser.add_argument("--output_prefix", 
                        default = "data/processed/fake_data", 
                        type= str
                       )
    parser.add_argument("--sample_nrows", 
                        default = 1000, 
                        type= int
                       )
    parser.add_argument("--inject_missing", 
                        default = True, 
                        type= bool
                       )
    args = parser.parse_args()
    main(args)