# 🔥 Ethanol-Water Distillation Heat Transfer Model  

This MATLAB-based project models the **heat transfer and mass separation** during the distillation of ethanol-water mixtures. Using condensation and boiling correlations, it calculates surface temperatures, heat flux, and mass transfer rates for different ethanol concentrations under varying steam conditions.  

---

## 🎯 Objective  

To compute and analyse the **heat transfer dynamics** of an ethanol-water distillation column, focusing on:  

- Calculating **surface temperatures** on steam and mixture sides of the plate.  
- Determining **heat transfer coefficients** for condensation and boiling.  
- Computing **heat flux** across the system.  
- Estimating **steam condensation rates** and **mixture vaporisation rates**.  
- Validating results across multiple ethanol concentrations and steam temperatures.  

---

## 🛠️ Tech Stack  

- **MATLAB** (iterative computation + tabular data handling)  
- **Custom functions** (heat transfer coefficient, Rayleigh/Nusselt numbers, viscosity, etc.)  
- **Experimental reference data** (for ethanol-water mixtures & steam)  

---

## 🧠 Key Concepts  

- Heat transfer by **condensation** (Gerstmann & Griffith correlation).  
- Heat transfer by **nucleate boiling** (Kutateladze correlation).  
- Iterative solution of coupled non-linear equations.  
- Use of **thermophysical property tables** for ethanol-water mixtures and steam.  
- Linear interpolation for intermediate steam temperatures.  

---

## ⚙️ Methodology  

1. **Problem Setup**  
   - Steam at temperatures: **150 °C, 180 °C, and 200 °C + K-number offset**.  
   - Mixture concentrations: **15%, 60%, 80%, 95.6% ethanol by weight**.  
   - Plate: Stainless steel (1 mm thick, k = 16.7 W/m·K).  

  <img src="Screenshot 2025-09-05 135907.png" height="70%" width="70%" />

2. **Equations Applied**  
   - Newton’s Law of Cooling (steam → wall).  
   - Fourier’s Law (conduction through plate).  
   - Boiling correlation (wall → mixture).  



3. **Iterative Solution**  
   - Start with an initial guess for wall temperature on the steam side.  
   - Use ΔT_sub to calculate Rayleigh & Nusselt numbers → condensation h.  
   - Compute heat flux q″ (constant across system).  
   - Use q″ to calculate boiling-side h_mix and T_w,mix.  
   - Refine guess for T_w,water until convergence.  


---

## 📊 Results  

Final values were computed for each mixture & steam temperature.  

- **Outputs per run**:  
  - T_w,water (°C)  
  - T_w,mix (°C)  
  - h_water (W/m²·K)  
  - h_mix (W/m²·K)  
  - q″ (W/m²)  
  - ṁ_water (kg/m²·hr)  
  - ṁ_mix (kg/m²·hr)  

<img src="Screenshot 2025-09-05 135957.png" height="70%" width="70%" />

---

## 🚧 Challenges  

- **Coupled equations**: Only 3 governing equations for 5 unknowns, requiring iterative interdependence.  
- **Calibration**: Heat transfer coefficients vary strongly with ethanol concentration.  
- **Runtime stability**: Iterations needed careful tolerance to converge reliably.  
- **Interpolation**: Steam properties at non-tabulated temperatures required linear interpolation.  

---

## ✅ Outcomes  

- Successfully computed **surface temperatures, heat flux, and transfer coefficients**.  
- Demonstrated ethanol concentration significantly affects boiling-side h.  
- Showed how **steam temperature increases improve distillation efficiency**.  
- Built a reusable MATLAB framework for multicomponent distillation modelling.  

---

## 📁 Repository Structure  

