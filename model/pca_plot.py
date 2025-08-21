import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.font_manager as fm

if __name__ == '__main__':
    generation_df = pd.read_csv('data/generation_inhibitor_10000_property.csv')
    inhibitor_df = pd.read_csv('data/inhibitor2_property.csv')


    combined_df = pd.concat([generation_df, inhibitor_df])


    X = combined_df[['logp', 'mw', 'tpsa', 'hbd', 'hba', 'rotbonds']].values

    font_path = r'times.ttf'
    prop = fm.FontProperties(fname=font_path, size=18, weight='bold')
    legend_prop = fm.FontProperties(fname=font_path, size=16, weight='bold')
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)


    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    generated_color = '#1f77b4'  
    fine_tuned_color = '#ff7f0e'  

    plt.figure(figsize=(8, 6))
   
    plt.scatter(X_pca[:len(generation_df), 0], X_pca[:len(generation_df), 1], c=generated_color, label='generated', alpha=0.6, s=30)
    plt.scatter(X_pca[len(generation_df):, 0], X_pca[len(generation_df):, 1], c=fine_tuned_color, label='fine-tuned', alpha=0.6, s=30)
    
    plt.title('')
    plt.xlabel('Principal Component 1', fontproperties=prop)
    plt.ylabel('Principal Component 2', fontproperties=prop)
    
    plt.legend(prop=prop, frameon=False)
    
    plt.grid(False)
    plt.tight_layout() 
    plt.savefig('PCA_2.svg',bbox_inches='tight')
    plt.show()
