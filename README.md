# Cost-Effectiveness of RSV Vaccination in U.S. Adults Aged 75–84  
### A Decision-Analytic and Machine Learning–Enhanced Modeling Study

This repository contains a complete and reproducible cost-effectiveness analysis (CEA) of **respiratory syncytial virus (RSV) vaccination** in **U.S. adults aged 75–84 years**. The study combines a **Markov health economic model** with **machine learning–based risk calibration** to evaluate both **clinical impact** and **economic value** of RSV vaccination from a **U.S. healthcare payer perspective**.

📄 **Final project report:** `CEA.pdf`  
📦 **Programming language:** R  
🧠 **Methods:** Decision analysis + Machine Learning  
📊 **Model outputs:** ICER, QALYs, CEAC, PSA

---

## ✅ Key Findings

| Outcome | Result |
|----------|--------|
| Incremental cost per person | **–$45.20** |
| Incremental QALYs per person | **+0.0000593** |
| ICER | **Dominant** (cost-saving + health gain) |
| Population savings (100,000 people) | **$4.52 million** |
| PSA probability cost-effective | **63–65%** at standard WTP thresholds |

---

## 🔧 Model Overview

| Feature | Description |
|----------|-------------|
| Model type | Cohort Markov |
| Health states | Healthy (H), Outpatient RSV (O), Hospitalized RSV (I), Death (D) |
| Time horizon | 15 years |
| Perspective | U.S. healthcare payer |
| Discounting | 3% (costs & QALYs) |
| Vaccine type | Arexvy-like waning VE |
| Sensitivity | Probabilistic Sensitivity Analysis (PSA, 2000 draws) |
| ML component | Logistic regression risk calibration |

---

## 📂 Repository Structure

