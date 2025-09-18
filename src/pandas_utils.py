import pandas as pd

def assign_with_per_row_fn(df, fn, new_cols, *args, **kwargs):
    return df.assign(
        **pd.DataFrame(
        df.apply(fn, axis=1, args=tuple(*args), **kwargs).to_list(),
        columns=new_cols,
        index=df.index
        )
    )
