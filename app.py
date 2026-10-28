
import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw, Crippen

st.set_page_config(
    page_title="Mohamed Fahad | Chemical Library",
    page_icon="⚗️",
    layout="wide"
)

st.markdown("""
<style>
.stApp {
    background: linear-gradient(135deg, #080d18, #14253b);
    color: #f1f5f9;
}
[data-testid="stSidebar"] {
    background-color: #0b1525;
}
[data-testid="stMetric"] {
    background: #172a42;
    padding: 16px;
    border-radius: 12px;
    border: 1px solid #304760;
}
h1, h2, h3 {
    color: #7dd3fc !important;
}
.stCaption {
    color: #94a3b8;
}
</style>
""", unsafe_allow_html=True)

st.title("⚗️ Mohamed Fahad")
st.subheader("Cheminformatics Library")
st.caption("MOLECULAR EXPLORER • CHEMICAL PROPERTIES • STRUCTURES")
st.divider()

# Molecule name, SMILES, category, description
MOLECULES = [
    ("Water", "O", "Inorganic", "Universal solvent"),
    ("Hydrogen peroxide", "OO", "Inorganic", "Oxidizing compound"),
    ("Ammonia", "N", "Inorganic", "Nitrogen hydride"),
    ("Methane", "C", "Hydrocarbon", "Simplest alkane"),
    ("Ethane", "CC", "Hydrocarbon", "Two-carbon alkane"),
    ("Propane", "CCC", "Hydrocarbon", "Three-carbon alkane"),
    ("Butane", "CCCC", "Hydrocarbon", "Four-carbon alkane"),
    ("Ethene", "C=C", "Hydrocarbon", "Simple alkene"),
    ("Ethyne", "C#C", "Hydrocarbon", "Simple alkyne"),
    ("Cyclohexane", "C1CCCCC1", "Hydrocarbon", "Cyclic alkane"),
    ("Benzene", "c1ccccc1", "Aromatic", "Aromatic hydrocarbon"),
    ("Toluene", "Cc1ccccc1", "Aromatic", "Methylbenzene"),
    ("Naphthalene", "c1ccc2ccccc2c1", "Aromatic", "Fused aromatic rings"),
    ("Phenol", "Oc1ccccc1", "Phenol", "Aromatic hydroxy compound"),
    ("Aniline", "Nc1ccccc1", "Amine", "Aromatic amine"),
    ("Anisole", "COc1ccccc1", "Ether", "Aromatic ether"),
    ("Nitrobenzene", "O=[N+]([O-])c1ccccc1", "Aromatic", "Nitro compound"),
    ("Pyridine", "n1ccccc1", "Heterocycle", "Nitrogen heterocycle"),
    ("Pyrrole", "c1cc[nH]c1", "Heterocycle", "Five-membered ring"),
    ("Furan", "o1cccc1", "Heterocycle", "Oxygen heterocycle"),
    ("Thiophene", "c1ccsc1", "Heterocycle", "Sulfur heterocycle"),
    ("Methanol", "CO", "Alcohol", "One-carbon alcohol"),
    ("Ethanol", "CCO", "Alcohol", "Two-carbon alcohol"),
    ("Propanol", "CCCO", "Alcohol", "Three-carbon alcohol"),
    ("Ethylene glycol", "OCCO", "Alcohol", "Two hydroxyl groups"),
    ("Glycerol", "OCC(O)CO", "Alcohol", "Three hydroxyl groups"),
    ("Diethyl ether", "CCOCC", "Ether", "Organic solvent"),
    ("Dimethyl ether", "COC", "Ether", "Simple ether"),
    ("Acetone", "CC(=O)C", "Ketone", "Common solvent"),
    ("Formaldehyde", "C=O", "Aldehyde", "Simplest aldehyde"),
    ("Acetaldehyde", "CC=O", "Aldehyde", "Simple aldehyde"),
    ("Formic acid", "O=CO", "Carboxylic acid", "One-carbon acid"),
    ("Acetic acid", "CC(=O)O", "Carboxylic acid", "Two-carbon acid"),
    ("Lactic acid", "CC(O)C(=O)O", "Carboxylic acid", "Hydroxy acid"),
    ("Benzoic acid", "O=C(O)c1ccccc1", "Carboxylic acid", "Aromatic acid"),
    ("Citric acid", "OC(=O)CC(O)(CC(O)=O)C(O)=O",
     "Carboxylic acid", "Organic acid"),
    ("Ethyl acetate", "CCOC(C)=O", "Ester", "Ester solvent"),
    ("Urea", "NC(N)=O", "Amide", "Nitrogen compound"),
    ("Acetamide", "CC(N)=O", "Amide", "Simple amide"),
    ("Glycine", "NCC(=O)O", "Amino acid", "Simplest amino acid"),
    ("Alanine", "CC(N)C(=O)O", "Amino acid", "Simple amino acid"),
    ("Glucose", "C(C1C(C(C(C(O1)O)O)O)O)O",
     "Carbohydrate", "Simple sugar"),
    ("Fructose", "C(C(C(C(C(=O)CO)O)O)O)O",
     "Carbohydrate", "Simple sugar"),
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O",
     "Pharmaceutical", "Acetylsalicylic acid"),
    ("Paracetamol", "CC(=O)Nc1ccc(O)cc1",
     "Pharmaceutical", "Acetaminophen"),
    ("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
     "Pharmaceutical", "Anti-inflammatory medicine"),
    ("Caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
     "Pharmaceutical", "Common stimulant"),
    ("Nicot ine", "CN1CCC[C@H]1c1cccnc1",
     "Alkaloid", "Naturally occurring alkaloid"),
    ("Vitamin C", "O=C1OC(C(O)CO)C(O)=C1O",
     "Biomolecule", "Ascorbic acid"),
    ("Dopamine", "NCCc1ccc(O)c(O)c1",
     "Biomolecule", "Catecholamine"),
    ("Serine", "N[C@@H](CO)C(=O)O",
     "Amino acid", "Protein amino acid"),
    ("Tryptophan", "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",
     "Amino acid", "Aromatic amino acid"),
]

@st.cache_data
def make_database():
    rows = []

    for name, smiles, category, description in MOLECULES:
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
            "Fraction Csp3": round(
                rdMolDescriptors.CalcFractionCSP3(mol), 3
            ),
            "Category": category,
            "Description": description,
            "SMILES": Chem.MolToSmiles(mol),
        })

    return pd.DataFrame(rows)


df = make_database()

# Sidebar filters
st.sidebar.title("🔬 Molecule Search")

search = st.sidebar.text_input(
    "Search by name, formula or SMILES"
)

categories = ["All"] + sorted(df["Category"].unique().tolist())
selected_category = st.sidebar.selectbox("Category", categories)

weight_range = st.sidebar.slider(
    "Molecular weight (g/mol)",
    min_value=0.0,
    max_value=float(df["Molecular Weight"].max()),
    value=(0.0, float(df["Molecular Weight"].max()))
)

filtered = df[
    df["Name"].str.contains(search, case=False, na=False)
    | df["Formula"].str.contains(search, case=False, na=False)
    | df["SMILES"].str.contains(search, case=False, na=False)
]

filtered = filtered[
    filtered["Molecular Weight"].between(
        weight_range[0], weight_range[1]
    )
]

if selected_category != "All":
    filtered = filtered[
        filtered["Category"] == selected_category
    ]

# Summary
m1, m2, m3 = st.columns(3)
m1.metric("Total Molecules", len(df))
m2.metric("Search Results", len(filtered))
m3.metric("Categories", df["Category"].nunique())

st.divider()

tab1, tab2, tab3 = st.tabs([
    "🧪 Explore",
    "📊 Compare",
    "📥 Download"
])

with tab1:
    st.header("Molecular Database")

    st.dataframe(
        filtered[[
            "Name", "Formula", "Molecular Weight",
            "LogP", "TPSA", "Category"
        ]],
        use_container_width=True,
        hide_index=True
    )

    if not filtered.empty:
        chosen = st.selectbox(
            "Choose a molecule for full details",
            filtered["Name"].tolist()
        )

        item = df[df["Name"] == chosen].iloc[0]
        mol = Chem.MolFromSmiles(item["SMILES"])

        left, right = st.columns([1, 1.4])

        with left:
            st.subheader(item["Name"])
            st.image(
                Draw.MolToImage(mol, size=(450, 350)),
                caption="2D Molecular Structure"
            )

        with right:
            st.subheader(item["Formula"])
            st.write(item["Description"])
            st.write("**Category:**", item["Category"])
            st.code(item["SMILES"], language=None)

            properties = [
                "Molecular Weight",
                "Exact Mass",
                "LogP",
                "TPSA",
                "H-Bond Donors",
                "H-Bond Acceptors",
                "Rotatable Bonds",
                "Heavy Atoms",
                "Ring Count",
                "Aromatic Rings",
                "Fraction Csp3",
            ]

            st.dataframe(
                pd.DataFrame({
                    "Property": properties,
                    "Value": [item[p] for p in properties]
                }),
                use_container_width=True,
                hide_index=True
            )
    else:
        st.warning("No molecules found. Change your search.")

with tab2:
    st.header("Compare Molecules")

    compare_names = st.multiselect(
        "Select up to 5 molecules",
        df["Name"].tolist(),
        default=["Ethanol", "Phenol", "Benzene", "Aspirin"]
    )

    if compare_names:
        comparison = df[
            df["Name"].isin(compare_names)
        ].set_index("Name")

        cols = [
            "Formula", "Molecular Weight", "Exact Mass",
            "LogP", "TPSA", "H-Bond Donors",
            "H-Bond Acceptors", "Ring Count"
        ]

        st.dataframe(
            comparison[cols].T,
            use_container_width=True
        )

        st.subheader("Molecular Weight Comparison")
        st.bar_chart(comparison["Molecular Weight"])

with tab3:
    st.header("Download Chemical Data")

    csv_data = filtered.to_csv(index=False).encode("utf-8")

    st.download_button(
        label="⬇️ Download CSV",
        data=csv_data,
        file_name="mohamed_fahad_molecules.csv",
        mime="text/csv"
    )

    st.caption("The download contains the currently filtered molecules.")

st.divider()

st.markdown("""
<div style="text-align:center; padding:18px; color:#9CA3AF;">
    <h3 style="color:#E2E8F0;">Mohamed Fahad</h3>
    <p>Cheminformatics Library | Molecule Explorer</p>
</div>
""", unsafe_allow_html=True)
