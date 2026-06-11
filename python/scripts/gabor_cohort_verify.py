import glob, os, numpy as np, pandas as pd
OUT="/ceph/margrie/laura/goggle_gabor/out"
fs=sorted(glob.glob(os.path.join(OUT,"*.parquet")))
print(f"{len(fs)} parquet files")
df=pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)
DPP=111.6/400.0
df["sf_tok"]=(df["cloud"].str.extract(r"_sf(\d+p\d+)_")[0].str.replace("p",".").astype(float)/DPP).round(3)
df["th_tok"]=(np.degrees(df["cloud"].str.extract(r"theta(-?\d+p\d+)")[0].str.replace("p",".").astype(float))%180).round()
print(f"total rows {len(df):,} | clouds {df.cloud.nunique()} | RFs {df[['probe','cluster','rf_type']].drop_duplicates().shape[0]} | edge-RF rows {int((df.edge>0).sum()):,}")
print(f"frames/cloud/RF: {df.groupby(['cloud','probe','cluster','rf_type']).size().unique()}")
print(f"NaNs: sf {df.sf_cpd.isna().sum()} or {df.or_deg.isna().sum()}")
print("\nSF recovery by token (mean obs vs token):")
for tok,g in df.groupby("sf_tok"):
    print(f"  token {tok:.3f}: obs {g.sf_cpd.mean():.4f} ± {g.sf_cpd.std():.4f}  (n={len(g):,})")
def cmean(d): return np.degrees(np.angle(np.mean(np.exp(1j*2*np.radians(d))))/2)%180
print("\nOR circular-mean by theta token:")
for tok,g in df.groupby("th_tok"):
    print(f"  token {int(tok):3d}: circ-mean {cmean(g.or_deg.values):.1f}  (n={len(g):,})")
print(f"\nconcentration (reliability): median {df.concentration.median():.1f}, "
      f"RFs with conc<2: {df[df.concentration<2][['probe','cluster','rf_type']].drop_duplicates().shape[0]}")
