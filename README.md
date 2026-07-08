# 🌡️ Dynamique des Vagues de Chaleur dans le Bassin Versant de la Bénoué (1981–2024)
### *Caractérisation spatio-temporelle, classification et modélisation probabiliste des événements thermiques extrêmes au Nord-Cameroun*

[![Licence: MIT](https://img.shields.io/badge/Licence-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version R](https://img.shields.io/badge/R-%3E%3D4.3-blue.svg)](https://www.r-project.org/)
[![Données ERA5](https://img.shields.io/badge/Données-ERA5%20ECMWF-green.svg)](https://cds.climate.copernicus.eu/)
[![Statut](https://img.shields.io/badge/Statut-Chapitre%201%20Thèse%20%E2%80%94%20Finalisé-brightgreen.svg)]()
[![Langue](https://img.shields.io/badge/Langue-Français-blue.svg)]()

---

## 📌 Présentation générale

Ce dépôt documente l'ensemble de la chaîne analytique du **Chapitre 1** d'une thèse de
doctorat en climatologie-épidémiologie au Cameroun.

À partir de **44 années de données climatiques multi-sources validées** (1981–2024), ce
chapitre livre la première caractérisation spatio-temporelle complète des vagues de chaleur
dans le **bassin versant de la Bénoué** (73 232 km², Nord-Cameroun) — une région sahélienne
abritant 3,5 millions d'habitants avec une infrastructure sanitaire fortement limitée.

**Résultat central :** la fréquence des vagues de chaleur a augmenté de **+341 %** entre
1981–2002 et 2003–2024, sous l'influence dominante de l'**Oscillation Multidécennale de
l'Atlantique** (AMO, ρ = +0,581, R² = 25,4 %), l'ENSO s'avérant non significatif
(ρ = −0,037) sur les températures extrêmes après correction correcte des données.

---

## 🔬 Contributions scientifiques originales

- **Classification KOUPNA II (2026) :** premier système gradué à 3 niveaux de sévérité
  des vagues de chaleur (C1 Normale / C2 À Risque / C3 Dangereuse) calibré sur les seuils
  physiologiques sahéliens et l'épidémiologie obstétricale — une première pour le Sahel
  camerounais
- **Analyse GEV bayésienne** (Ratio-of-Uniforms, N = 50 000) : distribution de Fréchet
  confirmée (ξ = +0,174, probabilité postérieure = 91,7 %) — T₁₀ = 46,0°C, T₅₀ = 48,3°C
- **Chaînes de Markov bayésiennes** aux échelles annuelle et journalière avec prévisions
  probabilistes à 14 jours opérationnelles pour les systèmes d'alerte précoce sanitaire
- **Cartographie ERA5 sur 253 pixels** agrégée en score composite d'exposition (SH) par
  commune GADM — identification de Taïbong, Guider et Figuil comme hotspots prioritaires
- **Rupture structurelle** 2001–2004 confirmée indépendamment par les tests de Pettitt
  et Bai-Perron sur trois stations météorologiques
- **Hiérarchie téléconnexions :** AMO > IOD > ENSO — résultat contre-intuitif démontré
  empiriquement pour la première fois pour les températures extrêmes sahéliennes

---

## 📂 Structure du dépôt

```
benoue-vagues-chaleur-1981-2024/
├── data/
│   ├── brutes/           # Données GSOD et GHCNd — NOAA (non publiées)
│   ├── ERA5/             # Fichiers NetCDF — 253 pixels, grille 0,25°
│   └── traitees/         # vagues_avec_classification.rds (N = 38 événements)
├── scripts/
│   ├── Module_01_10/     # Nettoyage, validation triple source
│   ├── Module_11_19/     # Détection, classification, tendances Mann-Kendall
│   ├── Module_20/        # Téléconnexions (ENSO, AMO, IOD/DMI)
│   ├── Module_21/        # Validation externe (3 stations)
│   └── Module_22/        # Cartographie spatiale (ERA5 × HydroSHEDS × GADM)
├── demos/
│   ├── demo_gev_bayesien.R           # Script démo GEV bayésien
│   ├── demo_mann_kendall_markov.R    # Script démo MK + Markov bayésien
│   └── README_demos.md               # Instructions d'exécution
├── outputs/
│   ├── figures/          # Figures publication (PNG 300 dpi)
│   └── tableaux/         # Exports CSV de tous les résultats
└── README.md
```

---

## 🛠️ Environnement technique (Tech Stack)

| Catégorie | Outils et packages |
|-----------|-------------------|
| **Langage principal** | R ≥ 4.3 |
| **Données climatiques** | ERA5-ECMWF (`ncdf4`, `terra`), GSOD-NOAA, GHCNd-NOAA |
| **Analyse spatiale** | `sf`, `terra`, HydroSHEDS Niveau 5, GADM v4.1 |
| **Tests de tendance et rupture** | `Kendall`, `trend` — Mann-Kendall, Sen, Pettitt, Bai-Perron |
| **Valeurs extrêmes** | `revdbayes`, `extRemes` — GEV bayésien, Ratio-of-Uniforms |
| **Modélisation Markov** | `MCMCpack` — matrices de transition bayésiennes (Dirichlet) |
| **Visualisation** | `ggplot2`, `ggpubr`, `ggspatial`, `patchwork` |
| **Manipulation des données** | `tidyverse`, `lubridate`, `stringr` |
| **Rapports et tableaux** | R Markdown, `flextable`, `officer` |
| **🤖 Assistance IA** | Les scripts R et le code analytique ont été développés avec l'assistance de l'**Intelligence Artificielle (outil LLM)**. L'ensemble des résultats numériques a été vérifié systématiquement par l'auteur sur les fichiers de données brutes avant intégration dans la thèse. |

---

## 📊 Résultats clés en un coup d'œil

```
Bassin          : Bénoué, Nord-Cameroun (73 232 km²)
Période         : 1981–2024 (44 ans, 16 071 jours-station)
Vagues P95      : N = 38  (seuil P95 = 38,5°C | OMM 2015)
Classification  : C1 = 31 (81,6 %) | C2 = 6 (15,8 %) | C3 = 1 (2,6 %)
Fréquence       : 0,32/an (1981–2002) → 1,41/an (2003–2024) = +341 %
Rupture         : 2001–2004  (Pettitt p < 0,01 + Bai-Perron confirmé)
Tendances MK    : τ = +0,452*** (HWF) | +0,600*** (J35) | +0,572*** (NT25)
GEV bayésien    : μ = 43,62°C | σ = 0,896 | ξ = +0,174 (Fréchet, 91,7 %)
Niveaux retour  : T10 = 46,0°C | T20 = 46,9°C | T50 = 48,3°C [IC95 : 46,5–56,0°C]
AMO             : ρ = +0,581*** | β = 4,61 vagues/°C | R² = 25,4 %
ENSO            : ρ = −0,037 (ns) — non significatif après correction des anomalies
Spatial ERA5    : τ moyen = +0,493 sur 253 pixels | pente Sen = +0,271°C/décennie
Vague C3 (1998) : 87,5 % de couverture du bassin | anomalie spatiale moyenne = +6,39°C
Record absolu   : 49,5°C — 19 février 2003 (Saison Sèche Froide, station GSOD Garoua)
```

---

## 🔍 Scripts de démonstration reproductibles

Deux scripts autonomes sont disponibles dans `/demos/` pour permettre
aux recruteurs, évaluateurs et pairs de vérifier les méthodes analytiques.

| Script | Méthodes illustrées |
|--------|---------------------|
| `demo_gev_bayesien.R` | GEV bayésien · Ratio-of-Uniforms · prior MDI · Niveaux de retour |
| `demo_mann_kendall_markov.R` | Mann-Kendall · Pente de Sen · Pettitt · Markov bayésien |

> Les deux scripts utilisent des **données simulées uniquement** et ne nécessitent
> que les packages R standards. Temps d'exécution < 2 minutes.

```r
# Installation des packages nécessaires
install.packages(c("revdbayes", "ggplot2", "patchwork",
                   "Kendall", "trend", "MCMCpack",
                   "tidyverse", "lubridate"))

# Exécution des scripts de démonstration
source("demos/demo_gev_bayesien.R")
source("demos/demo_mann_kendall_markov.R")
```

---

## 📖 Citation

```bibtex
@phdthesis{KoupnaII2026,
  author  = {Koupna II, Higelin Saint-Clair},
  title   = {Dynamique spatio-temporelle des vagues de chaleur et complications
             hypertensives de la grossesse (pré-éclampsie et éclampsie) :
             cartographie des facteurs de risque et étude des stratégies
             d'adaptation locales dans le bassin versant de la Bénoué
             (Nord-Cameroun)},
  year    = {2026},
  note    = {Chapitre 1 : Caractérisation spatio-temporelle des vagues de chaleur
             dans le bassin versant de la Bénoué, 1981--2024}
}
```

---

## 🔗 Sources de données

| Source | Utilisation | Accès |
|--------|-------------|-------|
| GSOD-NOAA | Données station primaires | https://www.ncdc.noaa.gov/cdo-web/ |
| GHCNd-NOAA | Contrôle qualité | https://www.ncdc.noaa.gov/ghcnd-data-access |
| ERA5-ECMWF | Analyse spatiale 253 pixels | https://cds.climate.copernicus.eu/ |
| HydroSHEDS N5 | Délimitation bassin versant | https://www.hydrosheds.org/ |
| GADM v4.1 | Limites administratives 39 communes | https://gadm.org/ |

---

## 👤 Auteur

**KOUPNA II Higelin Saint-Clair**
Doctorant indépendant — Climatologie-Épidémiologie
Chef de Service des Diplômes, Curricula et Recherche, Université de Garoua (FSEG)
Nord-Cameroun

---

## 🏥 Portée et valorisation

Ce dépôt s'inscrit dans une thèse sur les interactions climat-santé au Sahel camerounais.
Le **Chapitre 2** abordera les liens épidémiologiques entre la dynamique des vagues de
chaleur documentée ici et les complications hypertensives de la grossesse (pré-éclampsie
et éclampsie) dans le bassin versant de la Bénoué, à partir des données du HRG de Garoua
(2000–2024).

**Applications opérationnelles directes :**
- Système d'alerte précoce chaleur à 14 jours pour les districts sanitaires
- Cartographie de vulnérabilité thermique par commune pour les interventions ONG
- Révision du Plan National d'Adaptation aux Changements Climatiques (PNACC) Cameroun
- Outil KKRAIA+ (en développement) — interface décisionnelle climat-santé

---

*Ce projet est en libre accès. Les données brutes de station ne sont pas publiées
en raison de droits NOAA. Les fichiers ERA5 sont accessibles via le portail
Copernicus Climate Data Store (inscription gratuite requise).*

---

---

## 📬 Contact & Collaborations

I am available for **consulting missions, expert assessments and research
partnerships** in climate-health, early warning systems and climate
adaptation for Sub-Saharan Africa.

**Areas of expertise:**
- Extreme heat characterization and classification (Sahelian context)
- Climate-maternal health linkages (pre-eclampsia, eclampsia)
- Operational vulnerability mapping for health districts and NGOs
- Probabilistic forecasting and climate early warning systems

📩 **higeo.saintclair@gmail.com / gabigeo2@yahoo.fr**
🔗 LinkedIn: www.linkedin.com/in/koupna-ii-higelin-saint-clair-182b7922
