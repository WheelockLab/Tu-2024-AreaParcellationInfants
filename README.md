This is the repository for paper in submission: "[Early Life Neuroimaging: The Generalizability of Cortical Area Parcellations Across Development](https://www.biorxiv.org/content/10.1101/2024.09.09.612056v1.full)"
**Contact: Jiaxin Cindy Tu (tu.j@wustl.edu)**

**Abstract**

<p>The cerebral cortex comprises discrete cortical areas that form during development. Accurate area parcellation in neuroimaging studies enhances statistical power and comparability across studies. The formation of cortical areas is influenced by intrinsic embryonic patterning as well as extrinsic inputs, particularly through postnatal exposure. Given the substantial changes in brain volume, microstructure, and functional connectivity during the first years of life, we hypothesized that cortical areas in 1-to-3-year-olds would exhibit major differences from those in neonates and progressively resemble adults as development progresses.
Here, we parcellated the cerebral cortex into putative areas using local functional connectivity gradients in 92 toddlers at 2 years old. We demonstrated high reproducibility of these cortical regions across 1-to-3-year-olds in two independent datasets. The area boundaries in 1-to-3-year-olds were more similar to adults than neonates. While the age-specific group parcellation fitted better to the underlying functional connectivity in individuals during the first 3 years, adult area parcellations might still have some utility in developmental studies, especially in children older than 6 years. Additionally, we provided connectivity-based community assignments of the parcels, showing fragmented anterior and posterior components based on the strongest connectivity, yet alignment with adult systems when weaker connectivity was included. </p>

**File Organization**

**./Code** contains the code to generate and evaluate parcel performances.<br>
**./Adult and Infant Parcellations in 32k-fsLR** are the parcellations listed in Figure 3A and Table 2.<br>
**./Final 2-yr Parcellation** (Tu 326) is the area parcellation and boundary map for 92 healthy term 2-year-olds in the eLABE study in Figure 3A.<br>
**./Finer age-binned parcels** is the area parcellations generated with smaller age windows using the Baby Connectome Project (BCP) data in Figure 4A.<br>
**./NetworkCommunityAssignments** is the community assignment labels of the Tu 326 areas in Figure 6.<br>
**./2-yr Parcellation Different Resolutions** are the parcellations with different merging thresholds to generate different number of area parcels in Supplementary Figure 3.<br>
