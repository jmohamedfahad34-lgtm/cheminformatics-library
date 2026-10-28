
import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw
from rdkit.Chem import Crippen
from io import BytesIO

st.set_page_config(
    page_title="MoleculeVerse | Chemical Library",
    page_icon="⚗️",
    layout="wide"
)

st.markdown("""
<style>
.stApp {
    background: linear-gradient(135deg, #071426, #10233e);
    color: #f1f5f9;
}
[data-testid="stMetric"] {
    background: #142842;
    padding: 18px;
    border-radius: 12px;
    border: 1px solid #27435e;
}
section[data-testid="stSidebar"] {
    background: #0b1729;
}
h1, h2, h3 { color: #7dd3fc; }
</style>
""", unsafe_allow_html=True)

st.title("⚗️ MoleculeVerse")
st.caption("CHEMICAL EXPLORER • STRUCTURES • MOLECULAR DESCRIPTORS")

# Name, SMILES, category, short description
DATA = [
    ("Water", "O", "Inorganic", "Universal solvent"),
    ("Ethanol", "CCO", "Alcohol", "Alcohol used in solvents"),
    ("Methanol", "CO", "Alcohol", "Simple alcohol"),
    ("Propanol", "CCCO", "Alcohol", "Three-carbon alcohol"),
    ("Acetone", "CC(=O)C", "Ketone", "Common organic solvent"),
    ("Benzene", "c1ccccc1", "Aromatic", "Aromatic hydrocarbon"),
    ("Toluene", "Cc1ccccc1", "Aromatic", "Aromatic solvent"),
    ("Phenol", "Oc1ccccc1", "Phenol", "Aromatic hydroxy compound"),
    ("Aniline", "Nc1ccccc1", "Amine", "Aromatic amine"),
    ("Anisole", "COc1ccccc1", "Ether", "Aromatic ether"),
    ("Pyridine", "n1ccccc1", "Heterocycle", "Nitrogen heterocycle"),
    ("Pyrrole", "c1cc[nH]c1", "Heterocycle", "Five-membered heterocycle"),
    ("Furan", "o1cccc1", "Heterocycle", "Oxygen heterocycle"),
    ("Thiophene", "c1ccsc1", "Heterocycle", "Sulfur heterocycle"),
    ("Formic acid", "C(=O)O", "Carboxylic acid", "Simple carboxylic acid"),
    ("Acetic acid", "CC(=O)O", "Carboxylic acid", "Vinegar acid"),
    ("Benzoic acid", "O=C(O)c1ccccc1", "Carboxylic acid", "Aromatic acid"),
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O", "Pharmaceutical", "Acetylsalicylic acid"),
    ("Paracetamol", "CC(=O)Nc1ccc(O)cc1", "Pharmaceutical", "Acetaminophen"),
    ("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "Pharmaceutical", "Anti-inflammatory medicine"),
    ("Caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O", "Pharmaceutical", "Common stimulant"),
    ("Urea", "NC(=O)N", "Amide", "Nitrogen-containing compound"),
    ("Glycine", "NCC(=O)O", "Amino acid", "Simplest amino acid"),
    ("Glucose", "C(C1C(C(C(C(O1)O)O)O)O)O", "Carbohydrate", "Simple sugar"),
    ("Fructose", "C(C(C(C(C(=O)CO)O)O)O)O", "Carbohydrate", "Simple sugar"),
    ("Lactic acid", "CC(O)C(=O)O", "Carboxylic acid", "Hydroxy acid"),
    ("Citric acid", "C(C(=O)O)C(CC(=O)O)(C(=O)O)O", "Carboxylic acid", "Organic acid"),
    ("Naphthalene", "c1ccc2ccccc2c1", "Aromatic", "Fused aromatic rings"),
    ("Cyclohexane", "C1CCCCC1", "Hydrocarbon", "Saturated cyclic hydrocarbon"),
    ("Ethene", "C=C", "Hydrocarbon", "Simple alkene"),
    ("Ethyl acetate", "CCOC(=O)C", "Ester", "Common ester solvent"),
    ("Formaldehyde", "C=O", "Aldehyde", "Simplest aldehyde"),
    ("Acetaldehyde", "CC=O", "Aldehyde", "Simple aldehyde"),
    ("Dimethyl ether", "COC", "Ether", "Simple ether"),
    ("Nitrobenzene", "O=[N+]([O-])c1ccccc1", "Aromatic", "Nitro aromatic compound"),
]

@st.cache_data
def build_database():
    rows = []
    for name, smiles, category, description in DATA:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        rows.append({
            "Name": name,
            "Formula": rdMolDescriptors.CalcMolFormula(mol),
            "Molecular Weight": round(Descriptors.MolWt(mol), 3),
            "Exact Mass": round(Descriptors.ExactMolWt(mol), 4),
            "LogP": round(Crippen.MolLogP(mol), 3),
            "TPSA": round(Descriptors.TPSA(mol), 3),
            "H-Bond Donors": rdMolDescriptors.CalcNumHBD(mol),
            "H-Bond Acceptors": rdMolDescriptors.CalcNumHBA(mol),
            "Rotatable Bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "Heavy Atoms": mol.GetNumHeavyAtoms(),
            "Ring Count": rdMolDescriptors.CalcNumRings(mol),
            "Aromatic Rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "Fraction Csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
            "Category": category,
            "Description": description,
            "SMILES": Chem.MolToSmiles(mol),
        })
    return pd.DataFrame(rows)

df = build_database()

st.sidebar.title("🔬 Explorer")
query = st.sidebar.text_input("Search molecule")
categories = ["All"] + sorted(df["Category"].unique().tolist())
category = st.sidebar.selectbox("Category", categories)

max_mw = float(df["Molecular Weight"].max())
mw_range = st.sidebar.slider(
    "Molecular weight (g/mol)",
    0.0, max_mw, (0.0, max_mw)
)

filtered = df[
    df["Name"].str.contains(query, case=False, na=False)
    | df["Formula"].str.contains(query, case=False, na=False)
    | df["SMILES"].str.contains(query, case=False, na=False)
]

filtered = filtered[
    filtered["Molecular Weight"].between(*mw_range)
]

if category != "All":
    filtered = filtered[filtered["Category"] == category]

a, b, c = st.columns(3)
a.metric("Molecules", len(df))
b.metric("Search results", len(filtered))
c.metric("Categories", df["Category"].nunique())

tab1, tab2, tab3 = st.tabs([
    "🧪 Molecule Explorer",
    "📊 Compare",
    "📥 Export"
])

with tab1:
    st.subheader("Molecular Database")
    st.dataframe(
        filtered[[
            "Name", "Formula", "Molecular Weight",
            "LogP", "TPSA", "Category"
        ]],
        use_container_width=True,
        hide_index=True
    )

    if not filtered.empty:
        selected = st.selectbox(
            "Select a molecule",
            filtered["Name"].tolist()
        )
        row = df[df["Name"] == selected].iloc[0]
        mol = Chem.MolFromSmiles(row["SMILES"])

        left, right = st.columns([1, 1.5])

        with left:
            st.subheader(row["Name"])
            st.image(
                Draw.MolToImage(mol, size=(420, 320)),
                caption="2D molecular structure"
            )

        with right:
            st.markdown(f"### {row['Formula']}")
            st.write(row["Description"])
            st.write("**Category:**", row["Category"])
            st.write("**SMILES:**", row["SMILES"])

            props = [
                "Molecular Weight", "Exact Mass", "LogP",
                "TPSA", "H-Bond Donors", "H-Bond Acceptors",
                "Rotatable Bonds", "Heavy Atoms", "Ring Count",
                "Aromatic Rings", "Fraction Csp3"
            ]
            st.dataframe(
                pd.DataFrame({
                    "Property": props,
                    "Value": [row[p] for p in props]
                }),
                hide_index=True,
                use_container_width=True
            )
    else:
        st.info("No molecules found. Try another search.")

with tab2:
    st.subheader("Compare Molecules")
    names = st.multiselect(
        "Choose up to 5 molecules",
        df["Name"].tolist(),
        default=["Ethanol", "Phenol", "Benzene", "Aspirin"]
    )
    if names:
        comparison = df[df["Name"].isin(names)].set_index("Name")
        cols = [
            "Formula", "Molecular Weight", "LogP",
            "TPSA", "H-Bond Donors", "H-Bond Acceptors",
            "Rotatable Bonds", "Ring Count"
        ]
        st.dataframe(
            comparison[cols].T,
            use_container_width=True
        )
        st.bar_chart(comparison["Molecular Weight"])

with tab3:
    st.subheader("Download your results")
    csv = filtered.to_csv(index=False).encode("utf-8")
    st.download_button(
        "⬇️ Download filtered CSV",
        data=csv,
        file_name="molecule_library.csv",
        mime="text/csv"
    )
    st.caption(
        "Descriptors are computed from the SMILES in this dataset. "
        "Values are computational estimates, not experimental measurements."
    )

st.divider()
st.caption(
    "MoleculeVerse • Educational Cheminformatics Library | "
    "Created by Mohamed Fahad J, Joshipa I & Saivarshini"
)
