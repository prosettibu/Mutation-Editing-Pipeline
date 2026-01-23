import pandas as pd
import requests
import time
import json

def check_pathogenic(mutation):
    try:
        if mutation.get('rsid') and mutation['rsid'] != '' and str(mutation['rsid']) != 'nan':
            rsid = str(mutation['rsid']).replace('rs', '')
            
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            search_params = {
                'db': 'clinvar',
                'term': f'rs{rsid}[RS]',
                'retmode': 'json'
            }
            
            response = requests.get(search_url, params=search_params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                
                if 'esearchresult' in data and data['esearchresult'].get('idlist'):
                    clinvar_ids = data['esearchresult']['idlist']
                    
                    for clin_id in clinvar_ids[:3]:
                        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                        summary_params = {
                            'db': 'clinvar',
                            'id': clin_id,
                            'retmode': 'json'
                        }
                        
                        time.sleep(0.4)
                        response = requests.get(summary_url, params=summary_params, timeout=10)
                        
                        if response.status_code == 200:
                            summary = response.json()
                            
                            if 'result' in summary and clin_id in summary['result']:
                                variant_data = summary['result'][clin_id]
                                
                                clin_sig = None
                                if 'clinical_significance' in variant_data:
                                    clin_sig_data = variant_data['clinical_significance']
                                    if isinstance(clin_sig_data, dict):
                                        clin_sig = clin_sig_data.get('description', '')
                                    elif isinstance(clin_sig_data, str):
                                        clin_sig = clin_sig_data
                                
                                if not clin_sig and 'germline_classification' in variant_data:
                                    germline = variant_data['germline_classification']
                                    if isinstance(germline, dict):
                                        clin_sig = germline.get('description', '')
                                    elif isinstance(germline, str):
                                        clin_sig = germline
                                
                                if clin_sig:
                                    sig_lower = clin_sig.lower()
                                    
                                    if 'pathogenic' in sig_lower and 'benign' not in sig_lower:
                                        return True, 0.9
                                    elif 'benign' in sig_lower:
                                        return False, 0.1
                                    elif 'uncertain' in sig_lower or 'conflict' in sig_lower:
                                        return None, 0.5
        
        chrom = str(mutation['chrom']).replace('chr', '')
        
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        search_params = {
            'db': 'clinvar',
            'term': f'{mutation["gene"]}[gene] AND {chrom}[chr] AND {mutation["position"]}[chrpos37]',
            'retmode': 'json',
            'retmax': 10
        }
        
        time.sleep(0.4)
        response = requests.get(search_url, params=search_params, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            
            if 'esearchresult' in data and data['esearchresult'].get('idlist'):
                clinvar_ids = data['esearchresult']['idlist']
                
                for clin_id in clinvar_ids[:3]:
                    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                    summary_params = {
                        'db': 'clinvar',
                        'id': clin_id,
                        'retmode': 'json'
                    }
                    
                    time.sleep(0.4)
                    response = requests.get(summary_url, params=summary_params, timeout=10)
                    
                    if response.status_code == 200:
                        summary = response.json()
                        
                        if 'result' in summary and clin_id in summary['result']:
                            variant_data = summary['result'][clin_id]
                            
                            clin_sig = None
                            
                            if 'clinical_significance' in variant_data:
                                clin_sig_data = variant_data['clinical_significance']
                                if isinstance(clin_sig_data, dict):
                                    clin_sig = clin_sig_data.get('description', '')
                                else:
                                    clin_sig = str(clin_sig_data)
                            
                            if not clin_sig and 'germline_classification' in variant_data:
                                germline = variant_data['germline_classification']
                                if isinstance(germline, dict):
                                    clin_sig = germline.get('description', '')
                                else:
                                    clin_sig = str(germline)
                            
                            if clin_sig:
                                sig_lower = clin_sig.lower()
                                
                                if 'pathogenic' in sig_lower and 'benign' not in sig_lower:
                                    return True, 0.9
                                elif 'benign' in sig_lower:
                                    return False, 0.1
                                elif 'uncertain' in sig_lower or 'conflict' in sig_lower:
                                    return None, 0.5
        
        return False, 0.2
        
    except Exception:
        return None, 0.0

def process_csv(filename):
    df = pd.read_csv(filename)
    results = []
    
    for index, row in df.iterrows():
        mutation = {
            'gene': row['gene'],
            'chrom': row['chromosome'],
            'old_letter': row['ref'],
            'new_letter': row['alt'],
            'position': row['position'],
            'rsid': row.get('rsid', None)
        }
        
        is_bad, score = check_pathogenic(mutation)
        
        results.append({
            'gene': mutation['gene'],
            'position': f"chr{mutation['chrom']}:{mutation['position']}",
            'mutation': f"{mutation['old_letter']}→{mutation['new_letter']}",
            'rsid': mutation.get('rsid', ''),
            'pathogenic': is_bad,
            'score': score
        })
        
        if index < len(df) - 1:
            time.sleep(0.5)
    
    return pd.DataFrame(results)

if __name__ == "__main__":
    results = process_csv('mutations.csv')
    print(results.to_string(index=False))
    print(f"\nPathogenic: {sum(results['pathogenic'] == True)}")
    print(f"Benign: {sum(results['pathogenic'] == False)}")
    print(f"Uncertain: {sum(results['pathogenic'].isna())}")
    results.to_csv('pathogenicity_results.csv', index=False)