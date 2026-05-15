"""
	This script builds on the perl version by A.Roldan.

"""

import pathlib
import sympy as sp
import numpy as np
from scipy.optimize import least_squares
from Symbols_def import temp, constants
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
from Kinetics import RConstants





'''



SEPARATE SPECIES NAMES LETERS AND NUMBERS TO MAKE THE NUMBERS SUBSCRIPT IN PLOTS



GENERATE DATA AND PLOTS OF ir PROPORTIONAL TO THE CONCENTRATIONS


'''










class Diagnostics:
	def __init__(self, systems, processes, species, equation_factors, surfaces, adsorbates_by_surface,
							gas_number, surface_expr, experiment_data):
		temp_list = [float(i) for i in set(experiment_data[1:, 0].tolist())]
		time_list = [float(i) for i in set(experiment_data[1:, 1].tolist())]
		n_surfaces = len(surface_expr.keys())
		concentrations = {}
		for row in experiment_data[1:, :]:
			te = str(row[0])
			ti = str(row[1])
			values = row[2:-n_surfaces]     # removes the surface covarages
			if te not in concentrations:
				concentrations[te] = {}
			concentrations[te][ti] = [float(i) for i in values]

		functions, derivatives = Diagnostics.get_functions(systems, processes, species, equation_factors, surface_expr)
		clipping = (surfaces, adsorbates_by_surface, gas_number)

		''' IR evolution '''
		#Diagnostics.ir_evolution(systems, temp_list, time_list, concentrations)
		''' Rates '''
		Diagnostics.rates(temp_list, time_list, concentrations, functions)
		''' Degree of Equilibrium Control: thermodynamics '''
		Diagnostics.equilibrium_control(temp_list, time_list, concentrations, clipping, functions, derivatives)
		''' Degree of Rate/Selectivity Control: kinetics '''
		Diagnostics.assess_drc_validity(processes, temp_list, time_list, concentrations, clipping, functions,
										derivatives)

	@staticmethod
	def get_functions(systems, processes, species, equation_factors, surface_expr):
		# --- Symbols ---
		species_syms = {s: sp.Symbol(s) for s in species}
		k_list = list(processes.keys())
		k_syms = {pr: sp.Symbol(f"k_{pr}") for pr in processes}
		k_sym_list = [k_syms[k] for k in k_list]
		k_index = {k: i for i, k in enumerate(k_list)}		# Map k index
		concentrations_vec = sp.Matrix([species_syms[s] for s in species])
		all_variables = (temp, *concentrations_vec, *k_sym_list)

		# --- Reaction Step pairing ---
		step_pairs = [(k_list[i], k_list[i + 1]) for i in range(0, len(k_list), 2)]
		if not step_pairs:
			raise ValueError("No forward/backward step pairs found")
		step_labels = []
		for pr_f, pr_b in step_pairs:
			step_labels.append(f"r_{{{pr_f}}} \\leftrightharpoons r_{{{pr_b}}}")

		# --- Build RHS ---
		f_list = []
		for name in species:
			dydt = 0
			for j, pr in enumerate(k_list):
				expr = k_syms[pr]
				for i, r in enumerate(processes[pr]['reactants']):
					if r in species_syms:
						expr *= species_syms[r] ** processes[pr]['rstoichio'][i]
					else:
						expr *= surface_expr[r] ** processes[pr]['rstoichio'][i]
				dydt += equation_factors[name][j] * expr
			f_list.append(dydt)
		f_vec = sp.Matrix(f_list)
		f_func = sp.lambdify(all_variables, f_vec, "numpy")

		# --- Products stoichiometry ---
		products = [name for name in systems if systems[name]['kind'] == 'molecule' and systems[name]['pressure0'] == 0]
		stoichi_pairs = np.array([[equation_factors[p][k_index[pr_f]] - equation_factors[p][k_index[pr_b]]
								   for (pr_f, pr_b) in step_pairs] for p in products])
		s_plus_pairs = np.array([[max(equation_factors[p][k_list.index(pr_f)], 0.0) for (pr_f, _) in step_pairs]
								 for p in products])

		# --- Rate constants ---
		k_func = {k: sp.lambdify(temp, processes[k]['krate0'].subs(constants), 'numpy') for k in k_list}

		# --- Forward & backward rates ---
		r_forward_exprs = []
		r_backward_exprs = []
		for pr_f, pr_b in step_pairs:
			expr_f = k_syms[pr_f]		# forward
			for i, r in enumerate(processes[pr_f]['reactants']):
				if r in species_syms:
					expr_f *= species_syms[r] ** processes[pr_f]['rstoichio'][i]
				else:
					expr_f *= surface_expr[r] ** processes[pr_f]['rstoichio'][i]
			expr_b = k_syms[pr_b]		# backward
			for i, r in enumerate(processes[pr_b]['reactants']):
				if r in species_syms:
					expr_b *= species_syms[r] ** processes[pr_b]['rstoichio'][i]
				else:
					expr_b *= surface_expr[r] ** processes[pr_b]['rstoichio'][i]
			r_forward_exprs.append(expr_f)
			r_backward_exprs.append(expr_b)
		r_forward_vec = sp.Matrix(r_forward_exprs)
		r_backward_vec = sp.Matrix(r_backward_exprs)
		r_forward_func = sp.lambdify(all_variables, r_forward_vec, "numpy")
		r_backward_func = sp.lambdify(all_variables, r_backward_vec, "numpy")
		r_func = (r_forward_func, r_backward_func)

		functions = (k_list, k_func, f_vec, f_func, r_func, step_pairs, products, stoichi_pairs, s_plus_pairs, step_labels)

		# --- Stable net rate (symbolic) ---
		beta_vec = sp.Matrix([r_backward_exprs[i] / r_forward_exprs[i] for i in range(len(step_pairs))])
		r_net_vec = sp.Matrix([r_forward_exprs[i] * (1 - beta_vec[i]) for i in range(len(step_pairs))])

		# --- Jacobians and derivates ---
		j_sym = f_vec.jacobian(concentrations_vec)
		j_func = sp.lambdify(all_variables, j_sym, "numpy")
		df_dk_sym = f_vec.jacobian(k_sym_list)
		df_dk_func = sp.lambdify(all_variables, df_dk_sym, "numpy")
		dr_net_dk_sym = r_net_vec.jacobian(k_sym_list)
		dr_net_dtheta_sym = r_net_vec.jacobian(concentrations_vec)
		dr_net_dk_func = sp.lambdify(all_variables, dr_net_dk_sym, "numpy")
		dr_net_dtheta_func = sp.lambdify(all_variables, dr_net_dtheta_sym, "numpy")
		dr_forward_dk_sym = r_forward_vec.jacobian(k_sym_list)
		dr_forward_dtheta_sym = r_forward_vec.jacobian(concentrations_vec)
		dr_forward_dk_func = sp.lambdify(all_variables, dr_forward_dk_sym, "numpy")
		dr_forward_dtheta_func = sp.lambdify(all_variables, dr_forward_dtheta_sym, "numpy")

		derivatives = (j_func, df_dk_func, dr_net_dk_func, dr_net_dtheta_func, dr_forward_dk_func, dr_forward_dtheta_func)

		return functions, derivatives

	@staticmethod
	def clip_species(surfaces_s, adsorbates_by_surface, gas_number, y):
		for s_sym in surfaces_s:
			s_name = str(s_sym)
			for idx in adsorbates_by_surface[s_name]:  # idx is the position of the adsorbate in species
				y[idx] = np.clip(y[idx], 0.0, 1.0)  # adsorbates between 0 and 1
		for idx in gas_number:
			y[idx] = np.clip(y[idx], 0.0, None)  # ensures gases > 0
		return y

	@staticmethod
	def rates(temp_list, time_list, concentrations, functions):
		if len(time_list) < 2:  # Defensive conditions
			raise RuntimeError("ODE solver returned too few time points in Experiments.rates.")
		dt = time_list[-1] - time_list[0]
		if dt == 0 or not np.isfinite(dt):
			raise RuntimeError("Invalid time interval in ODE solution in Experiments.rates.")

		k_list, k_func, f_vec, f_func, r_func, step_pairs, products, stoichi_pairs, s_plus_pairs, step_labels = functions
		r_forward_func, r_backward_func = r_func

		header = ["Temperature"] + [f'{pr_f}/{pr_b}' for pr_f, pr_b in step_pairs]
		data_ss = [header]
		data_avg = [header]
		for temp_num in [temp_list[0], temp_list[-1]]:
			k_values = np.array([k_func[k](temp_num) for k in k_list])
			r_forward = []
			r_backward = []
			for i in time_list:
				conc_t = np.array(concentrations[str(temp_num)][str(i)])
				all_vals_t = (temp_num, *conc_t, *k_values)

				r_forward.append(r_forward_func(*all_vals_t))
				r_backward.append(r_backward_func(*all_vals_t))
			r_forward = np.array(r_forward)
			r_backward = np.array(r_backward)

			# --- Evaluate net rate ---
			beta = r_backward / (r_forward + 1e-30)		# stable net
			r_net = r_forward * (1 - beta)

			if np.any(~np.isfinite(r_net)):
				raise RuntimeError("Net rate contains NaN or inf -- Diagnostics.rates.")

			# --- Rate at steady State ---
			n_tail = max(1, int(0.1 * len(r_net)))  # if rate_time !=0, get 10% of the last points
			ss = np.mean(r_net[-n_tail:], axis=0).ravel()
			data_ss.append([temp_num, *ss])

			# --- Average rate ---
			avg = np.trapz(r_net, time_list, axis=0) / dt  # Numpy versions > 2.0 uses "trapezoid" instead of "trapz"
			data_avg.append([temp_num, *avg.ravel()])

		Diagnostics.printdata(['./KINETICS/ISOTHERMAL/', "SteadyState_Net_Rates"], data_ss)    # temperature x processes
		Diagnostics.printdata(['./KINETICS/ISOTHERMAL/', "Average_Net_Rates"], data_avg)  # temperature x processes
		ylabel = "TOF ($molecule \\cdot site^{-1} \\cdot s^{-1}$)"
		RConstants.log_barplot(['./KINETICS/ISOTHERMAL/',"SteadyState Net Rates"], "Reaction Step", ylabel, data_ss,
							   step_labels, 0.5)
		RConstants.log_barplot(['./KINETICS/ISOTHERMAL/',"Average Net Rates"], "Reaction Step", ylabel, data_avg,
								step_labels, 0.5)

	@staticmethod
	def equilibrium_control(temp_list, time_list, concentrations, clipping, functions, derivatives, eps=1e-30):
		''' It answers: the overall rate change upon a shift of the equilibrium of step i:
		 Compute Degree of Equilibrium Control (DEC)
				 DEC_i = DRC_forward_i − DRC_backward_i
			DEC ≈ 0    | step not controlling equilibrium
			DEC >> 0   | pushing equilibrium forward increases rate
			DEC << 0   | pushing backward increases rate
			DEC dominant | thermodynamic bottleneck	'''
		k_list, k_func, f_vec, f_func, r_func, step_pairs, products, stoichi_pairs, s_plus_pairs, step_labels = functions
		r_forward_func, r_backward_func = r_func
		k_index = {k: i for i, k in enumerate(k_list)}		# Map k index
		j_func, df_dk_func, dr_net_dk_func, dr_net_dtheta_func, dr_forward_dk_func, dr_forward_dtheta_func = derivatives

		def residual(x):
			return np.array(f_func(temp_num, *x, *k_values), dtype=float).ravel()

		data_dec = np.array(["Temperature", *[f'{pr_f}/{pr_b}' for pr_f, pr_b in step_pairs]])
		for temp_num in [temp_list[0], temp_list[-1]]:
			conc_guess = np.array(concentrations[str(temp_num)][str(time_list[-1])], dtype=float).ravel()
			k_values = np.array([k_func[k](temp_num) for k in k_list], dtype=float)
			sol = least_squares(residual, conc_guess, bounds=(0, np.inf))
			#sol = root(residual, conc_guess, method='lm')
			conc = Diagnostics.clip_species(*clipping, sol.x)
			all_vals = (temp_num, *conc, *k_values)

			# --- Evaluate rates ---
			r_forward = np.array(r_forward_func(*all_vals), dtype=float).flatten()  # (Nk,)
			r_backward = np.array(r_backward_func(*all_vals), dtype=float).flatten()
			beta = r_backward / (r_forward + 1e-30)		# stable net
			r_net = r_forward * (1 - beta)
			r_vec = stoichi_pairs @ r_net	# Products rates

			# --- Evaluate derivatives ---
			jacobian = np.array(j_func(*all_vals), dtype=float)
			df_dk = np.array(df_dk_func(*all_vals), dtype=float)
			dtheta_dk = -np.linalg.pinv(jacobian) @ df_dk
			dr_net_dk = np.array(dr_net_dk_func(*all_vals), dtype=float)
			dr_net_dtheta = np.array(dr_net_dtheta_func(*all_vals), dtype=float)
			dr_forward_dk = np.array(dr_forward_dk_func(*all_vals), dtype=float)
			dr_forward_dtheta = np.array(dr_forward_dtheta_func(*all_vals), dtype=float)

			# --- Total derivatives ---
			dr_net_total = (stoichi_pairs @ dr_net_dk) + (stoichi_pairs @ dr_net_dtheta) @ dtheta_dk
			dr_forward_total = (s_plus_pairs @ dr_forward_dk) + (s_plus_pairs @ dr_forward_dtheta) @ dtheta_dk

			dr_dk_total = dr_net_total

			dec = np.zeros((r_vec.shape[0], len(step_pairs)))
			for j, (pr_f, pr_b) in enumerate(step_pairs):
				kf = k_values[k_index[pr_f]]
				kb = k_values[k_index[pr_b]]

				# --- DRC contributions ---
				safe_r = np.where(np.abs(r_vec) < eps, np.nan, r_vec)
				drc_f = (kf * dr_dk_total[:, k_index[pr_f]]) / safe_r
				drc_b = (kb * dr_dk_total[:, k_index[pr_b]]) / safe_r

				# --- DEC ---
				dec[:, j] = drc_f - drc_b
			data_dec = np.vstack((data_dec, np.column_stack([np.tile(temp_num, (len(dec), 1)), dec])))

		Diagnostics.printdata(['./KINETICS/ISOTHERMAL/', "Equilibrium_Control"], data_dec)  # temperature x processes
		Diagnostics.barplot(['./KINETICS/ISOTHERMAL/', "Equilibrium Control"], "Reaction Step", "DEC", data_dec,
		                 step_labels, 0.5)

	@staticmethod
	def assess_drc_validity(processes, temp_list, time_list, concentrations, clipping, functions, derivatives,
							eq_threshold=1e-3, flux_threshold=1e-6,	rc_threshold=0.1, eps=1e-30):
		'''	(1) classifies steps
			(2) evaluates whether the DRC is meaningful
			(3) identifies rate-controlling steps
		for DRC, should be in steady-state but NOT in thermodynamic equilibrium		'''
		k_list, k_func, f_vec, f_func, r_func, step_pairs, products, stoichi_pairs, s_plus_pairs, step_labels = functions
		r_forward_func, r_backward_func = r_func

		def residual(x):
			return np.array(f_func(temp_num, *x, *k_values), dtype=float).ravel()

		# --- Loop for Temperatures ---
		data_eq = [["Temperature"] + [f"{r1}/{r2}" for r1, r2 in step_pairs]]
		do_drc = []
		for temp_num in [temp_list[0], temp_list[-1]]:
			conc_guess = np.array(concentrations[str(temp_num)][str(time_list[-1])], dtype=float).ravel()
			k_values = np.array([k_func[k](temp_num) for k in k_list], dtype=float)
			sol = least_squares(residual, conc_guess, bounds=(0, np.inf))
			conc = Diagnostics.clip_species(*clipping, sol.x)
			all_vals = (temp_num, *conc, *k_values)

			# --- Evaluate rates ---
			r_forward = np.array(r_forward_func(*all_vals), dtype=float).flatten()  # (Nk,)
			r_backward = np.array(r_backward_func(*all_vals), dtype=float).flatten()
			beta = r_backward / (r_forward + 1e-30)		# stable net
			r_net = r_forward * (1 - beta)

			# --- Metrics ---
			distance_eq = np.abs(r_net) / (np.abs(r_forward) + eps)
			flux_importance = np.abs(r_net) / (np.sum(np.abs(r_net)) + eps)
			data_eq.append(np.insert(distance_eq, 0, temp_num).tolist())

			# --- Classification ---
			output = open("KINETICS/ISOTHERMAL/DRC_Assessment.dat", "a+")
			output.write("#\n# Classification of steps to asses the validity of Degree of Rate Control\n")
			output.write(f"# Temperature = {temp_num} K\n")
			output.write("#   Process \t\tClassification Flux_importance distance_eq\n")
			labels = []
			for i in range(len(step_pairs)):
				if flux_importance[i] < flux_threshold:
					label = "dormant"
				elif distance_eq[i] < eq_threshold:
					label = "quasi_equilibrated"
				elif flux_importance[i] > rc_threshold:
					label = "rate_controlling"
				else:
					label = "intermediate"
				labels.append(label)
				output.write(f"{step_pairs[i]} \t{label} \t{flux_importance[i]} \t{distance_eq[i]}\n")

			# --- Global indicators ---
			global_distance = np.sum(np.abs(r_net)) / (np.sum(np.abs(r_forward)) + eps)
			n_rc = sum(l == "rate_controlling" for l in labels)
			do_drc.append(not (global_distance < 1e-6 or n_rc == 0))
			output.write(f"# Global_distance: {global_distance}\n# Is the DRC recommended: {do_drc[-1]}\n")
			output.close()

			# Plot: species in the reaction pathway
			Diagnostics.species_reaction_pathway(processes, step_pairs, r_net, temp_num, flux_threshold)

		# Plot: equilibrium proximity
		title = "DRC Assessment"
		xlabel = "Reaction step"
		ylabel = "Distance from equilibrium"

		RConstants.log_barplot(["KINETICS/ISOTHERMAL/", title], xlabel, ylabel, data_eq, step_labels,0.5)

		if any(do_drc):
			Diagnostics.kinetic_control(temp_list, time_list, clipping, functions, derivatives, concentrations)

	@staticmethod
	def species_reaction_pathway(processes, step_pairs, r_net, temp_num, flux_threshold):
		''' The Reaction Pathway plot shows reaction connectivity, not exact atom-by-atom transformation
		- Each node is one species appearing in your mechanism
		- It tells you which pathways dominate under the current conditions for products formation
				--> The thicker edges correspond to higher net flux (rate_forward / rate_backward)
		- Often the catalyst site node becomes a giant hub because many reactions consume/produce it.
		The plot is a flux-weighted reaction pathway map showing which species feed into which others and
		which catalytic pathways dominate under the chosen operating conditions. '''
		try:
			g = nx.DiGraph()
			for i, (pr_f, pr_b) in enumerate(step_pairs):
				flux = r_net[i]
				if abs(flux) < flux_threshold:
					continue

				pr = pr_f if flux >= 0 else pr_b   # Determine dominant direction

				reactants = processes[pr]['reactants']
				products = processes[pr]['products']

				for r in reactants:	    # Build edges reactants → products
					for p in products:
						if r == p:
							continue
						weight = abs(flux)
						if g.has_edge(r, p):
							g[r][p]['weight'] += weight
						else:
							g.add_edge(r, p, weight=weight)

			if g.number_of_edges() == 0:
				print("No significant fluxes → empty pathway graph")
				return

			pos = nx.spring_layout(g, seed=0)
			edges = list(g.edges())
			weights = [g[u][v]['weight'] for u, v in edges]
			max_w = max(weights)
			fig = plt.figure(figsize=(10, 6), clear=True)
			nx.draw(g, pos, with_labels=True)
			nx.draw_networkx_edges(g, pos, width=[(w / max_w) * 5 for w in weights])
			plt.savefig(f"KINETICS/ISOTHERMAL/Reaction_Pathway_{temp_num}.svg", dpi=300, orientation='landscape',
			            transparent=True)
			plt.close(fig)

			# --- Extract numerical pathway table
			output = open("KINETICS/ISOTHERMAL/Reaction_Pathway.dat", "a")
			output.write("#\n# Flux-weighted pathway: each edge r → p carries a numerical weight representing\n")
			output.write("#   cumulative net reaction flux from reactant to product\n\n")  # headings
			output.write(f"# Temperature = {temp_num} K\n#     Process_species \t\tweighted flux\n")

			pathway_data = [(u, v, d['weight']) for u, v, d in g.edges(data=True)]
			pathway_data.sort(key=lambda x: x[2], reverse=True)
			for u, v, w in pathway_data:    # Dominant reaction pathways
				output.write(f"{u:>10} → {v:<10} \t{w:.4e}\n")
			output.close()
		except ImportError:
			print("networkx not installed → skipping pathway plot")

	@staticmethod
	def kinetic_control(temp_list, time_list, clipping, functions, derivatives, concentrations):
		''' The Degree of Rate Control (DRC), introduced by C. T. Campbell (J. Catal. 204, 520, 2001), quantifies how
		   sensitive the overall reaction rate is to each elementary step’s rate constant. '''
		''' DRC must be computed using a reaction rate, not a species balance.
		Correct choices:	rate of product formation
							rate of reactant consumption
							sum of elementary fluxes forming product '''
		k_list, k_func, f_vec, f_func, r_func, step_pairs, products, stoichi_pairs, s_plus_pairs, step_labels = functions
		r_forward_func, r_backward_func = r_func
		j_func, df_dk_func, dr_net_dk_func, dr_net_dtheta_func, dr_forward_dk_func, dr_forward_dtheta_func = derivatives

		def residual(x):
			return np.array(f_func(temp_num, *x, *k_values), dtype=float).ravel()

		# --- Results ---
		data_heading = ["Temperature"] + [f"Step_{pr_f}/{pr_b}" for pr_f, pr_b in step_pairs]
		drc_data = {p: [] for p in products}
		dsc_data = {p: [] for p in products}

		# --- Loop for Temperatures ---
		for temp_num in [temp_list[0], temp_list[-1]]:
			conc_guess = np.array(concentrations[str(temp_num)][str(time_list[-1])], dtype=float).ravel()
			k_values = np.array([k_func[k](temp_num) for k in k_list], dtype=float)
			sol = least_squares(residual, conc_guess, bounds=(0, np.inf))
			conc = Diagnostics.clip_species(*clipping, sol.x)
			all_vals = (temp_num, *conc, *k_values)

			# --- Evaluate rates ---
			r_forward = np.array(r_forward_func(*all_vals), dtype=float).flatten()  # (Nk,)
			r_backward = np.array(r_backward_func(*all_vals), dtype=float).flatten()
			beta = r_backward / (r_forward + 1e-30)		# stable net
			r_net = r_forward * (1 - beta)
			r_vec = stoichi_pairs @ r_net	# Products rates

			# --- Evaluate derivatives ---
			jacobian = np.array(j_func(*all_vals), dtype=float)
			df_dk = np.array(df_dk_func(*all_vals), dtype=float)
			dtheta_dk = -np.linalg.pinv(jacobian) @ df_dk
			dr_net_dk = np.array(dr_net_dk_func(*all_vals), dtype=float)
			dr_net_dtheta = np.array(dr_net_dtheta_func(*all_vals), dtype=float)
			dr_forward_dk = np.array(dr_forward_dk_func(*all_vals), dtype=float)
			dr_forward_dtheta = np.array(dr_forward_dtheta_func(*all_vals), dtype=float)

			# --- Total derivatives ---
			dr_net_total = (stoichi_pairs @ dr_net_dk) + (stoichi_pairs @ dr_net_dtheta) @ dtheta_dk
			dr_forward_total = (s_plus_pairs @ dr_forward_dk) + (s_plus_pairs @ dr_forward_dtheta) @ dtheta_dk

			# --- Hybrid weighting adapting far or close to equilibrium ---
			eps = 1e-12
			alpha = np.abs(r_vec) / (np.abs(r_vec) + np.sum(r_forward) + eps)
			# alpha = np.clip(np.abs(r_vec) / (np.sum(np.abs(r_forward)) + eps), 0.0, 1.0)     # more aggressive
			dr_dk_total = alpha[:, None] * dr_net_total + (1 - alpha[:, None]) * dr_forward_total
			mask = np.abs(r_vec) > 1e-20      # Mask
			r_vec = r_vec[mask]
			dr_dk_total = dr_dk_total[mask]
			filtered_products = [p for i, p in enumerate(products) if mask[i]]

			# --- DRC ---
			x = dr_dk_total / (r_vec[:, None] + 1e-30)
			for i, p in enumerate(filtered_products):
				row = [temp_num]
				for j, (pr_f, _) in enumerate(step_pairs):
					row.append(float(x[i, j]))
				drc_data[p].append(row)

			# --- DSC ---
			r_tot = np.sum(r_vec)
			if abs(r_tot) < 1e-20:
				continue

			s_vec = r_vec / r_tot
			kdr = dr_dk_total
			kdr_sum = np.sum(kdr, axis=0)
			dsc = (kdr - s_vec[:, None] * kdr_sum) / (s_vec[:, None] * r_tot)
			for i, p in enumerate(filtered_products):
				row = [temp_num]
				for j, (pr_f, _) in enumerate(step_pairs):
					row.append(float(dsc[i, j]))
				dsc_data[p].append(row)

		for p in products:
			data_drc = np.vstack([data_heading, drc_data[p]])
			Diagnostics.printdata(['./KINETICS/ISOTHERMAL/', f"{p}_Degree_of_Rate_Control"], data_drc)
			ylabel = f"$DRC_{{{p}}}$"
			Diagnostics.barplot(["KINETICS/ISOTHERMAL/", f"{p}_DRC"], "Reaction Step", ylabel,
			                    data_drc, step_labels, 0.5)
			data_dsc = np.vstack([data_heading, dsc_data[p]])
			Diagnostics.printdata(['./KINETICS/ISOTHERMAL/', f"{p}_Degree_of_Selectivity_Control"], data_dsc)
			ylabel = f"$DSC_{{{p}}}$"
			Diagnostics.barplot(["KINETICS/ISOTHERMAL/", f"{p}_DSC"], "Reaction Step", ylabel,
			                    data_dsc, step_labels, 0.5)

	@staticmethod
	def barplot(experiment, x_label, y_label, data, xticklabels, bar_width):
		icolour = ["b", "r", "c", "g", "m", "y", "grey", "olive", "brown", "pink", "darkgreen", "seagreen", "khaki",
		   "teal"]
		ipatterns = ["///", "...", "xx", "**", "\\", "|", "--", "++", "oo", "OO"]
		gap = 0.7   # gap between group of columns, e.g. processes
		data = np.asarray(data)
		data = np.vstack([data[0], data[1:][data[1:, 0].argsort()]])  # sort remaining rows by Temperatures
		if isinstance(xticklabels, list):
			labels = xticklabels
		else:
			labels = data[0][1:].copy()
		x = np.arange(len(labels)) * (1 + gap)
		temps = data[[1, -1], 0]	# temperatures from first and last data rows (skipping header row at index 0)
		y = data[[1, -1], 1:len(labels) + 1].astype(float)	# values from first and last data rows (skip temperature)
		y_min = min([min(y[n]) for n in range(len(y))])
		y_max = max([max(y[n]) for n in range(len(y))])
		y_limit = [y_min-0.15*abs(y_min), y_max+0.15*y_max]

		fig, ax1 = plt.subplots(figsize=(10, 6), clear=True)
		ax1.set_ylim(bottom=np.min(y) * 0.5, top=np.max(y) * 5)  # Prevents from hitting zero during autoscaling
		for n in range(len(y)):    # processes | the first column is temperature
			offset = (n - len(temps)/2) * bar_width + bar_width/2
			ax1.bar(x + offset, y[n], log=False,
					width=bar_width, hatch=ipatterns[n], color=icolour[n], label=f"{temps[n]} K", alpha=0.7)

		x_limit = [ax1.get_xlim()[0], ax1.get_xlim()[1]]
		ax1.plot(x_limit, [0, 0], "k-", lw=1.5)
		ax1.set_xlim(x_limit)
		ax1.set_xlabel(x_label, fontsize=18)
		ax1.set_xticks(x)
		ax1.tick_params(axis='x', labelrotation=45, labelsize=14, labelright=True, pad=-2.5)
		ax1.set_xticklabels([f"${i}$" for i in labels], rotation=45, ha="right", va="top")
		for tick, label in zip(ax1.get_xticks(), ax1.get_xticklabels()):
			label.set_horizontalalignment('center')
			#label.set_x(tick + 0.1)
		ax1.set_ylim(y_limit)
		ax1.set_ylabel(y_label, fontsize=18)
		ax1.tick_params(axis='y', rotation=0, labelsize=16)

		leg_lines, leg_labels = ax1.get_legend_handles_labels()
		legend = ax1.legend(leg_lines[0:len(temps)], leg_labels[0:len(temps)], loc='best', fontsize=16)
		plt.title(experiment[-1])
		fig.tight_layout()
		plt.ion()
		plt.savefig(experiment[0] + "_".join(experiment[-1].split()) + ".svg", dpi=300, orientation='landscape',
					transparent=True)

	@staticmethod
	def printdata(experiment, data):
		maxlen = [max([len(f"{data[r][c]}") + 2 for r in range(len(data))]) for c in range(len(data[0]))]
		folder = experiment[0]
		outputfile = folder + "/" + str(experiment[1]) + ".dat"
		if not pathlib.Path(folder).exists():
			pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
		output = open(outputfile, "a")
		output.write("#")
		for i in range(len(data[0])):
			output.write(" {val:>{wid}s}".format(wid=maxlen[i], val=str(data[0][i])))  # headings
		output.write("\n")
		for row in data[1:]:
			for i in range(len(row)):
				a = np.abs(float(row[i]))
				output.write(" {val:>{wid}.3{c}}".format(wid=maxlen[i], val= a, c='f' if 1e-3 < a < 1e3 or a == 0. else
				'e'))
			output.write("\n")
		output.close()
