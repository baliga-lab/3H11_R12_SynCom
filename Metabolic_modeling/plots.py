import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


class PlotFit:

    def __init__(self, model, data, od_to_biomass_coeff, cobra):
        self.data = data
        self.model = model
        self.cobra = cobra
        self.solutions_bio = []
        self.solutions_atp = {}
        self.od_to_biomass_coeff = od_to_biomass_coeff
        self.intervals = len(self.data['time']) - 1

    def fit(self):
        for i in range(self.intervals):
            print(i, i + 1)
            _t0 = self.data['time'][i]
            _t1 = self.data['time'][i + 1]
            _t = _t1 - _t0
            _ac0 = self.data['Acetate'][i]
            _ac1 = self.data['Acetate'][i + 1]

            _od600_0 = self.data['OD600'][i]
            _od600_1 = self.data['OD600'][i + 1]

            _no3_0 = self.data['NO3'][i]
            _no3_1 = self.data['NO3'][i + 1]

            _ac_uptake = (_ac1 - _ac0) / _t
            _no3_uptake = (_no3_1 - _no3_0) / _t
            _od600_change = (_od600_1 - _od600_0) / _t
            _biomass_change = _od600_change * self.od_to_biomass_coeff

            if _ac_uptake < 0 and _no3_uptake < 0:
                print('time', _t0, _t1)
                print(f'Acetate {_ac0} - {_ac1} ({_ac_uptake}/time)')
                print(f'NO3 {_no3_0} - {_no3_1} ({_no3_uptake}/time)')
                print(f'OD600 {_od600_0} - {_od600_1} ({_od600_change}/time, Biomass {_biomass_change}/time) ')

                self.model.objective = 'bio1'
                self.model.reactions.EX_cpd00209_e0.lower_bound = _no3_uptake
                self.model.reactions.EX_cpd00029_e0.lower_bound = _ac_uptake
<<<<<<< HEAD
                self.model.reactions.bio1.lower_bound = 0
                self.model.reactions.bio1.upper_bound = _biomass_change
=======
                self.model.reactions.bio1.bounds = (0, _biomass_change)
>>>>>>> 29bdfec6c037c5e19430fe11d9a3281e12bfb64c

                sol_bio = self.cobra.flux_analysis.pfba(self.model)
                self.solutions_bio.append(sol_bio)

                _biomass_err = (_biomass_change - sol_bio.fluxes['bio1']) ** 2
                _ac_err = (_ac_uptake - sol_bio.fluxes['EX_cpd00029_e0']) ** 2
                _no3_err = (_no3_uptake - sol_bio.fluxes['EX_cpd00209_e0']) ** 2

                print('error biomass:', _biomass_err, 'AC', _ac_err, 'NO3', _no3_err)

                self.model.objective = 'ATPM_c0'
<<<<<<< HEAD
                self.model.reactions.bio1.lower_bound = _biomass_change
                self.model.reactions.bio1.upper_bound = _biomass_change
=======
                self.model.reactions.bio1.bounds = (_biomass_change, _biomass_change)
>>>>>>> 29bdfec6c037c5e19430fe11d9a3281e12bfb64c

                sol_atp = self.cobra.flux_analysis.pfba(self.model)
                self.solutions_atp[(_t0, _t1)] = sol_atp

                _biomass_err = (_biomass_change - sol_atp.fluxes['bio1']) ** 2
                _ac_err = (_ac_uptake - sol_atp.fluxes['EX_cpd00029_e0']) ** 2
                _no3_err = (_no3_uptake - sol_atp.fluxes['EX_cpd00209_e0']) ** 2

                print('error biomass:', _biomass_err, 'AC', _ac_err, 'NO3', _no3_err)
<<<<<<< HEAD
                print(sol_atp.fluxes['ATPM_c0'])
=======
                print(sol_atp.fluxes['ATPM_c0'], _biomass_change/sol_atp.fluxes['ATPM_c0'])
>>>>>>> 29bdfec6c037c5e19430fe11d9a3281e12bfb64c
                print()

    @staticmethod
    def generate_total_acc_data(solutions, capture, base_values=None):
        if base_values is None:
            base_values = {}
        df_res = {k: [] for k in capture}
        df_res['time'] = []
        arr_acc = {k: base_values.get(k, 0) for k in capture}

        for (_t0, _t1), sol in solutions.items():
            time = _t0
            _t = _t1 - _t0
            for time_step in range(_t):
                for alias, rxn_id in capture.items():
                    df_res[alias].append(arr_acc[alias])
                    arr_acc[alias] += sol.fluxes[rxn_id]
                df_res['time'].append(time)
                time += 1

        return df_res

    def build_plot_data(self, od_t0=0.01093, acetate_t0=18.0313603, no3_t0=9.14450462, capture=None):
        if capture is None:
            capture = {
                'biomass': 'bio1',
                'acetate': 'EX_cpd00029_e0',
                'no3': 'EX_cpd00209_e0',
                'no2': 'EX_cpd00075_e0',
                'n2o': 'EX_cpd00659_e0',
                'no': 'EX_cpd00418_e0',
                'n2': 'EX_cpd00528_e0',
            }
        base_values = {
            'biomass': od_t0 * self.od_to_biomass_coeff,
            'acetate': acetate_t0,
            'no3': no3_t0,
        }
        _data_pred = PlotFit.generate_total_acc_data(self.solutions_atp, capture, base_values)

        return _data_pred

    def get_exp_data(self):
        _data_exp = {k: list(v.values()) for k, v in self.data.items()}
        _data_exp['biomass'] = [x * self.od_to_biomass_coeff for x in _data_exp['OD600']]
        _data_exp.keys()
        return _data_exp

    def plot(self, _data_pred, _data_exp):
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt

        color_acetate = 'black'
        color_no3 = 'blue'
        color_no2 = 'green'
        color_n2o = 'purple'
        color_no = 'orange'
        color_n2 = 'cyan'

        fig, ax = plt.subplots()
        ax2 = ax.twinx()
        fig.set_size_inches(25.7, 8.27)

        sns.lineplot(data=_data_pred, x='time', y='biomass', ax=ax2, color='sienna')
        sns.scatterplot(data=_data_exp, x='time', y='biomass', ax=ax2, color='sienna')

        sns.lineplot(data=_data_pred, x='time', y='acetate', ax=ax, color=color_acetate)
        sns.scatterplot(data=_data_exp, x='time', y='Acetate', ax=ax, color=color_acetate)
        sns.lineplot(data=_data_pred, x='time', y='no3', ax=ax, color=color_no3)
        sns.scatterplot(data=_data_exp, x='time', y='NO3', ax=ax, color=color_no3)
        sns.lineplot(data=_data_pred, x='time', y='no2', ax=ax, color=color_no2)
        sns.scatterplot(data=_data_exp, x='time', y='NO2', ax=ax, color=color_no2)
        sns.lineplot(data=_data_pred, x='time', y='n2o', ax=ax, color=color_n2o)
        sns.scatterplot(data=_data_exp, x='time', y='N2O', ax=ax, color=color_n2o)

        sns.lineplot(data=_data_pred, x='time', y='no', ax=ax, color=color_no)
        sns.lineplot(data=_data_pred, x='time', y='n2', ax=ax, color=color_n2)

        plt.legend(handles=[
            mpatches.Patch(color=color_acetate, label='Acetate'),
            mpatches.Patch(color=color_no3, label='NO3'),
            mpatches.Patch(color=color_no2, label='NO2'),
            mpatches.Patch(color=color_no, label='NO'),
            mpatches.Patch(color=color_n2o, label='N2O'),
            mpatches.Patch(color=color_n2, label='N2'),
            mpatches.Patch(label='Dot experimental values', facecolor=None, color='white'),
            mpatches.Patch(edgecolor='sienna', label='Biomass', facecolor='white'),
        ])
        ax.set_title('Concentration')
        ax.set_xlabel('Time')
        ax.set_ylabel('mM')
        ax2.set_ylabel('gDW')
        sns.despine(fig, ax)
<<<<<<< HEAD

=======
                
>>>>>>> 29bdfec6c037c5e19430fe11d9a3281e12bfb64c

class CommPlots:

    def __init__(self, model):
        self.linestyle_3h11 = ':'
        self.linestyle_r12 = '-.'
        self.color_acetate = 'black'
        self.model = model

    def generate_solutions(self, cobra, rates, expected_growth_60R_40A):
        """
        generates growth solutions and max ATPM of R12 and 3H11
        :param cobra:
        :param rates:
        :param expected_growth_60R_40A:
        :return:
        """
        self.model.medium = {'EX_cpd00067_e0': 100.0,
                             'EX_cpd00058_e0': 100.0,
                             'EX_cpd00001_e0': 1000.0,
                             'EX_cpd00971_e0': 34.9324073,
                             'EX_cpd00013_e0': 1000.0,
                             'EX_cpd00244_e0': 100.0,
                             'EX_cpd00205_e0': 1.3415688,
                             'EX_cpd00009_e0': 100.0,
                             'EX_cpd11574_e0': 100.0,
                             'EX_cpd00305_e0': 100.0,
                             'EX_cpd00048_e0': 100.0,
                             'EX_cpd00209_e0': 12.0,
                             'EX_cpd00254_e0': 100.0,
                             'EX_cpd03424_e0': 100.0,
                             'EX_cpd15574_e0': 100.0,
                             'EX_cpd10515_e0': 100.0,
                             'EX_cpd00149_e0': 100.0,
                             'EX_cpd00029_e0': 20.0,
                             'EX_cpd00034_e0': 100.0,
                             'EX_cpd00063_e0': 100.0,
                             'EX_cpd00030_e0': 100.0,
                             'EX_cpd00099_e0': 100.0}
        solution_array = {}

        self.model.objective_direction = 'max'
        for i in range(len(rates['EX_cpd00029_e0'])):
            self.model.objective = 'bio1'

            uptake_ac = rates['EX_cpd00029_e0'][i]
            uptake_no3 = rates['EX_cpd00209_e0'][i] if rates['EX_cpd00209_e0'][i] < 0 else 0
            uptake_no2 = rates['EX_cpd00075_e0'][i] if rates['EX_cpd00075_e0'][i] < 0 else 0
            # uptake_ac = rates['EX_cpd00029_e0'][i]
            _medium_up = {
                'EX_cpd00029_e0': uptake_ac,
                'EX_cpd00209_e0': uptake_no3,
                'EX_cpd00075_e0': uptake_no2,
            }
            for ex_id in _medium_up.keys():
                rxn_ex = self.model.reactions.get_by_id(ex_id)
                rxn_ex.lower_bound = _medium_up[ex_id]
                # print(i, rxn_ex.id, rxn_ex.lower_bound)
            sol = cobra.flux_analysis.pfba(self.model)
            solution_array[i] = {'solution': sol}
            extra_growth = sol.fluxes['bio1'] - expected_growth_60R_40A[i]

            print(i, 'GROWTH 0 ATPM, exceess', sol.fluxes['bio1'], extra_growth)

            if extra_growth > 0:
                self.model.reactions.bio1.lower_bound = expected_growth_60R_40A[i]
                self.model.objective = 'ATPM_cA'
                sol_atpm = cobra.flux_analysis.pfba(self.model)
                atpm_max_A = sol_atpm.fluxes['ATPM_cA']
                solution_array[i]['max_atpm_a'] = atpm_max_A
                self.model.objective = 'ATPM_cR'
                sol_atpm = cobra.flux_analysis.pfba(self.model)
                atpm_max_R = sol_atpm.fluxes['ATPM_cR']
                solution_array[i]['max_atpm_r'] = atpm_max_A
                print(i, 'MAX isolate ATPM', atpm_max_A, atpm_max_R)
            else:
                solution_array[i]['max_atpm_a'] = 0
                solution_array[i]['max_atpm_r'] = 0
            # solution_array[i] = {'solution': sol}

            self.model.reactions.bio1.lower_bound = 0
            self.model.reactions.ATPM_cA.lower_bound = 0
            self.model.objective = 'bio1'

        return solution_array

    def generate_growth_gap_solutions(self, cobra, expected_growth_60R_40A, rates, solution_array):
        """
        Generate ATPM solutions with minimum growth rate from expected growth data
        :param cobra:
        :param expected_growth_60R_40A:
        :param rates:
        :param solution_array:
        :return:
        """
        solution_exp = {}

        self.model.reactions.ATPM_cA.lower_bound = 0
        self.model.reactions.ATPM_cR.lower_bound = 0
        self.model.reactions.bio1.lower_bound = 0
        self.model.optimize()

        for i in range(len(expected_growth_60R_40A)):
            uptake_acetate = rates['EX_cpd00029_e0'][i]
            uptake_nitrate = rates['EX_cpd00209_e0'][i] if rates['EX_cpd00209_e0'][i] < 0 else 0
            uptake_nitrite = rates['EX_cpd00075_e0'][i] if rates['EX_cpd00075_e0'][i] < 0 else 0

            self.model.reactions.EX_cpd00029_e0.lower_bound = uptake_acetate
            self.model.reactions.EX_cpd00209_e0.lower_bound = uptake_nitrate
            self.model.reactions.EX_cpd00075_e0.lower_bound = uptake_nitrite

            self.model.reactions.ATPM_cA.lower_bound = solution_array[i]['max_atpm_a'] * 0.4
            self.model.reactions.bio1.lower_bound = expected_growth_60R_40A[i]
            self.model.objective = 'ATPM_cR'

            print('Acetate:', uptake_acetate,
                  'NO3', uptake_nitrate,
                  'NO2', uptake_nitrite,
                  'min community biomass', expected_growth_60R_40A[i], 'mATP 3H11',
                  self.model.reactions.ATPM_cA.lower_bound)

            solution = cobra.flux_analysis.pfba(self.model)

            solution_exp[i] = solution

            print('Found solution mATP R12', solution_exp[i].fluxes['ATPM_cR'], 'solution status', solution.status)
            
        return solution_exp

    @staticmethod
    def constants():
        atpm_max_cA_isolate = 9.67570638881224
        atpm_max_cR_isolate = 23.801193202545065
        time_steps = [23, 19, 11, 18, 24, 24, 23]
        expected_growth_60R_40A = [3.98504E-06, 4.40109E-05, 3.79833E-05, 0.00001748, 8.79E-06, 0.000008511,
                                   -1.46191E-06]
        rates = {
            # acetate ac
            'EX_cpd00029_e0': [-0.00796337, -0.077158311, -0.126205391, -0.047720639, -0.044917204, -0.0163845,
                               -0.018624539],
            # nitrate NO3
            'EX_cpd00209_e0': [-0.015765185, -0.260613674, -0.336007637, 0.000367934, -0.00050591, 0.000160971,
                               -0.00033594],
            # nitrite NO2
            'EX_cpd00075_e0': [-0.000300664, 0.145752713, 0.30446552, -0.112949487, -0.100448798, -0.06937439,
                               0.000601328],
            # Nitrous oxide N2O
            'EX_cpd00659_e0': [0.002927536, 0.007122807, -0.018121213, -3.70367E-05, 2.77775E-05, -0.000138889, 0]
        }

        return atpm_max_cA_isolate, atpm_max_cR_isolate, time_steps, expected_growth_60R_40A, rates

    @staticmethod
    def get_exp_syncom():
        exp_data_symcom = {
            'time': [0, 23, 42, 53, 71, 95, 119, 142],
            'biomass': [0.00008208, 0.000173736, 0.001009944, 0.00142776, 0.0017424, 0.00195336, 0.002157624, 0.002124],
            'acetate': [18.0313603, 17.8482028, 16.3821949, 14.9939356, 14.1349641, 13.0569512, 12.6637232, 12.2353588],
            'no3': [9.14450462, 8.78190537, 3.83024557, 0.13416156, 0.14078438, 0.12864254, 0.13250585, 0.12477923],
            'no2': [0.09977915, 0.09286388, 2.86216542, 6.21128614, 4.17819538, 1.76742424, 0.10243888, 0.11626942],
            'n2o': [0, 0.06733333, 0.20266667, 0.00333333, 0.00266667, 0.00333333, 0, 0]
        }
        df_exp = {
            'i': [],
            'acetate': [],
            'biomass': [],
            'no3': [],
            'no2': [],
            'n2o': []
        }
        for i in range(7):
            df_exp['i'].append(exp_data_symcom['time'][i])
            df_exp['biomass'].append(exp_data_symcom['biomass'][i])
            df_exp['acetate'].append(exp_data_symcom['acetate'][i])
            df_exp['no3'].append(exp_data_symcom['no3'][i])
            df_exp['no2'].append(exp_data_symcom['no2'][i])
            df_exp['n2o'].append(exp_data_symcom['n2o'][i])

        return df_exp

    @staticmethod
    def exp_atpm_pred(model, od_coeff, d):
        model.reactions.ATPM_c0.lower_bound = 0
        res = {x: [] for x in ['total_time', 'gdw', 'ac', 'no3', 'no2', 'n2o', 'atpm', 'atpm_per_gdw_per_h']}
        for i in range(len(d) - 1):
            d_start = d[i]
            d_end = d[i + 1]
            time = d_end['time'] - d_start['time']
            dt_acetate = (d_end['Acetate'] - d_start['Acetate']) / time
            dt_OD = (d_end['OD600'] - d_start['OD600']) / time
            dt_no3 = (d_end['NO3'] - d_start['NO3']) / time
            dt_no2 = (d_end['NO2'] - d_start['NO2']) / time
            dt_n2o = (d_end['N2O'] - d_start['N2O']) / time
            res['total_time'].append(time)
            res['gdw'].append(dt_OD * od_coeff)
            res['ac'].append(dt_acetate)
            res['no3'].append(dt_no3)
            res['no2'].append(dt_no2)
            res['n2o'].append(dt_n2o)
            max_atpm = None
            if dt_acetate < 0 and dt_no3 < 0:
                model.objective = 'ATPM_c0'
                model.reactions.bio1.lower_bound = dt_OD * od_coeff
                model.reactions.EX_cpd00209_e0.lower_bound = dt_no3
                model.reactions.EX_cpd00029_e0.lower_bound = dt_acetate
                try:
                    sol_max_atpm = model.optimize()
                    if sol_max_atpm.status == 'optimal':
                        max_atpm = sol_max_atpm.fluxes['ATPM_c0']
                except Exception as e:
                    pass

            res['atpm'].append(max_atpm)
            res['atpm_per_gdw_per_h'].append(max_atpm / time if max_atpm else None)

        return pd.DataFrame(res)

    def generate_total_acc_data(self, time_steps, solution_exp):
        df_array = {
            'i': [],
            'gdw_3h11': [],
            'gdw_r12': [],
            'gdw_total': [],
            'acetate': [],
            'NO3': [],
            'NO2': [],
            'N2O': [],
            'N2': [],
            'NO': [],
        }
        acc_r_growth_acc = 0.00005472
        acc_a_growth_acc = 0.00002736
        acc_acetate = 18.0313603
        acc_no3 = 9.14450462
        acc_no2 = 0.09977915
        acc_no = 0
        acc_n2o = 0
        acc_n2 = 0
        time = 0
        for t_index in range(len(time_steps)):
            sol_exp = solution_exp[t_index]
            growth_3H11 = sol_exp.fluxes['bio1_A']
            growth_R12 = sol_exp.fluxes['bio1_R']
            up_acetate = sol_exp.fluxes['EX_cpd00029_e0']
            up_no3 = sol_exp.fluxes['EX_cpd00209_e0']
            up_no2 = sol_exp.fluxes['EX_cpd00075_e0']
            up_n2 = sol_exp.fluxes['EX_cpd00528_e0']
            up_n2o = sol_exp.fluxes['EX_cpd00659_e0']
            up_no = sol_exp.fluxes['EX_cpd00418_e0']

            iterations = time_steps[t_index]
            for it in range(iterations):
                df_array['i'].append(time)
                df_array['gdw_3h11'].append(acc_a_growth_acc)
                df_array['gdw_r12'].append(acc_r_growth_acc)
                df_array['gdw_total'].append(acc_r_growth_acc + acc_a_growth_acc)
                df_array['acetate'].append(acc_acetate)
                df_array['NO3'].append(acc_no3)
                df_array['NO2'].append(acc_no2)
                df_array['N2O'].append(acc_n2o)
                df_array['N2'].append(acc_n2)
                df_array['NO'].append(acc_no)
                acc_r_growth_acc += growth_R12
                acc_a_growth_acc += growth_3H11
                acc_no3 += up_no3
                acc_no2 += up_no2
                acc_n2o += up_n2o
                acc_n2 += up_n2
                acc_no += up_no
                acc_acetate += up_acetate
                time += 1

        print('Biomass 3H11', df_array['gdw_3h11'][-1],
              'Biomass R12', df_array['gdw_r12'][-1],
              'Biomass Total', df_array['gdw_3h11'][-1] + df_array['gdw_r12'][-1],
              'Time', df_array['i'][-1],
              'Acetate', df_array['acetate'][-1],
              'NO3', df_array['NO3'][-1],
              'NO2', df_array['NO2'][-1],
              'N2O', df_array['N2O'][-1],
              'N2', df_array['N2'][-1])

        return df_array

    def plot_total_acc(self, df_array, df_exp):


        fig, ax = plt.subplots()
        ax2 = ax.twinx()
        fig.set_size_inches(25.7, 8.27)
        # sns.lineplot(data = df_array[0], x='i', y='value', hue='reaction', marker='o', sort = False, ax=ax, linestyle='-')
        sns.lineplot(data=df_array, x='i', y='gdw_3h11', ax=ax2, linestyle=self.linestyle_3h11, color='sienna')
        sns.lineplot(data=df_array, x='i', y='gdw_r12', ax=ax2, linestyle=self.linestyle_r12, color='sienna')
        sns.lineplot(data=df_array, x='i', y='gdw_total', ax=ax2, linestyle='--', color='sienna')
        sns.scatterplot(data=df_exp, x='i', y='biomass', ax=ax2, color='sienna')
        sns.lineplot(data=df_array, x='i', y='acetate', ax=ax, color=self.color_acetate)
        sns.scatterplot(data=df_exp, x='i', y='acetate', ax=ax, color=self.color_acetate)
        sns.lineplot(data=df_array, x='i', y='NO3', ax=ax, color='blue')
        sns.scatterplot(data=df_exp, x='i', y='no3', ax=ax, color='blue')
        sns.lineplot(data=df_array, x='i', y='NO2', ax=ax, color='green')
        sns.scatterplot(data=df_exp, x='i', y='no2', ax=ax, color='green')
        sns.lineplot(data=df_array, x='i', y='N2O', ax=ax, color='purple')
        sns.scatterplot(data=df_exp, x='i', y='n2o', ax=ax, color='purple')
        sns.lineplot(data=df_array, x='i', y='NO', ax=ax, color='orange', linestyle='--')
        sns.lineplot(data=df_array, x='i', y='N2', ax=ax, color='cyan')

        plt.legend(handles=[
            mpatches.Patch(color=self.color_acetate, label='Acetate'),
            mpatches.Patch(color='blue', label='NO3'),
            mpatches.Patch(color='green', label='NO2'),
            mpatches.Patch(color='orange', label='NO'),
            mpatches.Patch(color='purple', label='N2O'),
            mpatches.Patch(color='cyan', label='N2'),
            mpatches.Patch(label='Dot experimental values', facecolor=None, color='white'),
            mpatches.Patch(edgecolor='sienna', label='Biomass 3H11', linestyle=self.linestyle_3h11, facecolor='white'),
            mpatches.Patch(edgecolor='sienna', label='Biomass R12', linestyle=self.linestyle_r12, facecolor='white'),
            mpatches.Patch(edgecolor='sienna', label='Biomass Total', linestyle='--', facecolor='white'),
        ])
        ax.set_title('Totals')
        ax.set_xlabel('time')
        ax.set_ylabel('mM')
        ax2.set_ylabel('gDW')
        sns.despine(fig, ax)

    def generate_total_uptake_data(self, time_steps, solution_exp):
        df_array = {
            'i': [],
            'gdw_3h11': [],
            'gdw_r12': [],
            'atpm_3h11': [],
            'atpm_r12': [],
            'acetate': [],
            'NO3': [],
            'NO2': [],
            'N2O': [],
            'N2': [],
            'NO': [],
        }
        time = 0
        for t_index in range(len(time_steps)):
            sol_exp = solution_exp[t_index]
            growth_3H11 = sol_exp.fluxes['bio1_A']
            growth_R12 = sol_exp.fluxes['bio1_R']
            atpm_3H11 = sol_exp.fluxes['ATPM_cA']
            atpm_R12 = sol_exp.fluxes['ATPM_cR']
            print(atpm_3H11, atpm_R12)
            up_acetate = sol_exp.fluxes['EX_cpd00029_e0']
            up_no3 = sol_exp.fluxes['EX_cpd00209_e0']
            up_no2 = sol_exp.fluxes['EX_cpd00075_e0']
            up_n2 = sol_exp.fluxes['EX_cpd00528_e0']
            up_n2o = sol_exp.fluxes['EX_cpd00659_e0']
            up_no = sol_exp.fluxes['EX_cpd00418_e0']

            iterations = time_steps[t_index]
            for it in range(iterations):
                df_array['i'].append(time)
                df_array['gdw_3h11'].append(growth_3H11)
                df_array['gdw_r12'].append(growth_R12)
                df_array['atpm_3h11'].append(atpm_3H11)
                df_array['atpm_r12'].append(atpm_R12)
                df_array['acetate'].append(up_acetate)
                df_array['NO3'].append(up_no3)
                df_array['NO2'].append(up_no2)
                df_array['N2O'].append(up_n2o)
                df_array['N2'].append(up_n2)
                df_array['NO'].append(up_no)
                time += 1

        return df_array

    def plot_total_uptake(self, df_array):
        fig, ax = plt.subplots(3, 1, height_ratios=[3, 1, 1])
        fig.set_size_inches(25.7, 8.27)

        sns.lineplot(data=df_array, x='i', y='acetate', ax=ax[0], color=self.color_acetate)
        sns.lineplot(data=df_array, x='i', y='NO3', ax=ax[0], color='blue')
        sns.lineplot(data=df_array, x='i', y='NO2', ax=ax[0], color='green')
        sns.lineplot(data=df_array, x='i', y='N2O', ax=ax[0], color='purple')
        sns.lineplot(data=df_array, x='i', y='NO', ax=ax[0], color='orange', linestyle='--')
        sns.lineplot(data=df_array, x='i', y='N2', ax=ax[0], color='cyan')

        sns.lineplot(data=df_array, x='i', y='gdw_3h11', ax=ax[1], linestyle=self.linestyle_3h11,
                     color='sienna')
        sns.lineplot(data=df_array, x='i', y='gdw_r12', ax=ax[1], linestyle=self.linestyle_r12,
                     color='sienna')
        sns.lineplot(data=df_array, x='i', y='atpm_3h11', ax=ax[2], linestyle=self.linestyle_3h11,
                     color='dimgray')
        sns.lineplot(data=df_array, x='i', y='atpm_r12', ax=ax[2], linestyle=self.linestyle_r12,
                     color='dimgray')
        ax[0].set_xlabel('')
        ax[1].set_xlabel('')
        ax[2].set_xlabel('time')
        ax[0].set_ylabel('mM/gDW/h')
        ax[1].set_ylabel('Biomass 1/h')
        ax[2].set_ylabel('mATP mM/gDW/h')
        sns.despine(fig, ax[0])
        sns.despine(fig, ax[1])
        sns.despine(fig, ax[2])

    def generate_organism_uptake_data(self, time_steps, solution_exp, monitor=None):
        if monitor is None:
            monitor = {
                'acetate': 'cpd00029_e0',
                'no3': 'cpd00209_e0',
                'no2': 'cpd00075_e0',
                'no': 'cpd00418_e0',
                'n2o': 'cpd00659_e0',
                'n2': 'cpd00528_e0',
                
                'leu': 'cpd00107_e0',
            }
        df_array = []
        
        for i in range(3):
            _data_array = {k: [] for k in monitor}
            _data_array['i'] = []
            _data_array['growth'] = []
            df_array.append(_data_array)
            
        time = 0
        for t_index in range(len(time_steps)):
            sol_exp = solution_exp[t_index]
            growth_3H11 = sol_exp.fluxes['bio1_A']
            growth_R12 = sol_exp.fluxes['bio1_R']
            #leu_total, leu_acido, leu_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00107_e0', sol_exp))
            #ac_total, ac_acido, ac_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00029_e0', sol_exp))
            #no3_total, no3_acido, no3_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00209_e0', sol_exp))
            #no2_total, no2_acido, no2_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00075_e0', sol_exp))
            #no_total, no_acido, no_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00418_e0', sol_exp))
            #n2o_total, n2o_acido, n2o_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00659_e0', sol_exp))
            #sn2_total, n2_acido, n2_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, 'cpd00528_e0', sol_exp))
            
            iterations = time_steps[t_index]
            
            for alias, cpd in monitor.items():
                cpd_total, cpd_acido, cpd_rhoda = self._cpd_acc_to_a_r_t(self._get_cpd_acc(self.model, cpd, sol_exp))
                for it in range(iterations):
                    df_array[0][alias].append(cpd_total)
                    df_array[1][alias].append(cpd_acido)
                    df_array[2][alias].append(cpd_rhoda)

            #print(ac_total, ac_acido, ac_rhoda)

            
            for it in range(iterations):
                df_array[0]['i'].append(time)
                df_array[1]['i'].append(time)
                df_array[2]['i'].append(time)
                df_array[1]['growth'].append(sol_exp.fluxes['bio1_A'])
                df_array[2]['growth'].append(sol_exp.fluxes['bio1_R'])
                time += 1
                """
                df_array[0]['acetate'].append(ac_total)
                df_array[1]['acetate'].append(ac_acido)
                df_array[2]['acetate'].append(ac_rhoda)

                df_array[0]['no3'].append(no3_total)
                df_array[1]['no3'].append(no3_acido)
                df_array[2]['no3'].append(no3_rhoda)

                df_array[0]['no2'].append(no2_total)
                df_array[1]['no2'].append(no2_acido)
                df_array[2]['no2'].append(no2_rhoda)

                df_array[0]['no'].append(no_total)
                df_array[1]['no'].append(no_acido)
                df_array[2]['no'].append(no_rhoda)

                df_array[0]['n2o'].append(n2o_total)
                df_array[1]['n2o'].append(n2o_acido)
                df_array[2]['n2o'].append(n2o_rhoda)

                df_array[0]['n2'].append(n2_total)
                df_array[1]['n2'].append(n2_acido)
                df_array[2]['n2'].append(n2_rhoda)

                df_array[0]['leu'].append(leu_total)
                df_array[1]['leu'].append(leu_acido)
                df_array[2]['leu'].append(leu_rhoda)
                """
                

        return df_array

    def plot_organism_uptake_data(self, df_array, title='Organism Uptake/Secretion rates'):
        fig, ax_arr = plt.subplots(3, 1, height_ratios=[3, 1, 1])
        fig.set_size_inches(25.7, 8.27)
        ax = ax_arr[0]
        ax2 = ax_arr[1]
        ax3 = ax_arr[2]

        sns.lineplot(data=df_array[1], x='i', y='no3', ax=ax3, color='blue', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='no3', ax=ax3, color='blue', linestyle=self.linestyle_r12)

        sns.lineplot(data=df_array[1], x='i', y='no2', ax=ax3, color='green', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='no2', ax=ax3, color='green', linestyle=self.linestyle_r12)

        sns.lineplot(data=df_array[1], x='i', y='leu', ax=ax2, color='green', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='leu', ax=ax2, color='green', linestyle=self.linestyle_r12)

        sns.lineplot(data=df_array[1], x='i', y='acetate', ax=ax, color='black', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='acetate', ax=ax, color='black', linestyle=self.linestyle_r12)

        sns.lineplot(data=df_array[1], x='i', y='no', ax=ax, color='orange', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='no', ax=ax, color='orange', linestyle=self.linestyle_r12)

        sns.lineplot(data=df_array[1], x='i', y='n2o', ax=ax, color='purple', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='n2o', ax=ax, color='purple', linestyle=self.linestyle_r12)

        sns.lineplot(data=df_array[1], x='i', y='n2', ax=ax, color='cyan', linestyle=self.linestyle_3h11)
        sns.lineplot(data=df_array[2], x='i', y='n2', ax=ax, color='cyan', linestyle=self.linestyle_r12)
        ax.set_title(title)

        ax.legend(handles=[
            mpatches.Patch(edgecolor=self.color_acetate, label='Acetate (3H11)', facecolor='white',
                           linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor=self.color_acetate, label='Acetate (R12)', facecolor='white',
                           linestyle=self.linestyle_r12),
            mpatches.Patch(edgecolor='orange', label='NO (3H11)', facecolor='white', linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor='orange', label='NO (R12)', facecolor='white', linestyle=self.linestyle_r12),
            mpatches.Patch(edgecolor='purple', label='N2O (3H11)', facecolor='white', linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor='purple', label='N2O (R12)', facecolor='white', linestyle=self.linestyle_r12),
            mpatches.Patch(edgecolor='cyan', label='N2 (3H11)', facecolor='white', linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor='cyan', label='N2 (R12)', facecolor='white', linestyle=self.linestyle_r12),
        ])
        ax2.legend(handles=[
            mpatches.Patch(edgecolor='green', label='Leucine (3H11)', facecolor='white', linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor='green', label='Leucine (R12)', facecolor='white', linestyle=self.linestyle_r12),
        ])
        ax3.legend(handles=[
            mpatches.Patch(edgecolor='blue', label='NO3 (3H11)', facecolor='white', linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor='blue', label='NO3 (R12)', facecolor='white', linestyle=self.linestyle_r12),
            mpatches.Patch(edgecolor='green', label='NO2 (3H11)', facecolor='white', linestyle=self.linestyle_3h11),
            mpatches.Patch(edgecolor='green', label='NO2 (R12)', facecolor='white', linestyle=self.linestyle_r12),
        ])

        ax.set_ylabel('mM/gDW/h')
        ax2.set_ylabel('mM/gDW/h')
        ax3.set_ylabel('mM/gDW/h')
        ax.set_xlabel('')
        ax2.set_xlabel('')
        ax3.set_xlabel('time')
        sns.despine(fig, ax)
        sns.despine(fig, ax2)
        sns.despine(fig, ax3)

    def generate_leucine_efflux_data(self, cobra, solution_atpm, expected_growth_60R_40A, rates, per_atpm_to_leucine):
        solution_leucine = {}
        for i in range(len(expected_growth_60R_40A)):
            solution_leucine[i] = None
            self.model.reactions.EX_cpd00029_e0.lower_bound = rates['EX_cpd00029_e0'][i]
            self.model.reactions.EX_cpd00209_e0.lower_bound = rates['EX_cpd00209_e0'][i] if \
            rates['EX_cpd00209_e0'][
                i] < 0 else 0
            self.model.reactions.EX_cpd00075_e0.lower_bound = rates['EX_cpd00075_e0'][i] if \
            rates['EX_cpd00075_e0'][
                i] < 0 else 0
            self.model.reactions.ATPM_cA.lower_bound = solution_atpm[i].fluxes['ATPM_cA'] * (1 - per_atpm_to_leucine)
            self.model.reactions.ATPM_cR.lower_bound = solution_atpm[i].fluxes['ATPM_cR'] * (1 - per_atpm_to_leucine)
            self.model.reactions.bio1.lower_bound = expected_growth_60R_40A[i]
            self.model.objective = 'LeuE_cA'
            self.model.objective_direction = 'min'
            try:
                solution_leucine[i] = cobra.flux_analysis.pfba(self.model)
                print(i, solution_leucine[i].fluxes['LeuE_cA'])
            except Exception as ex:
                pass
        self.model.objective_direction = 'max'
        self.model.objective = 'bio1'
        return solution_leucine

    @staticmethod
    def _get_cpd_acc(model, cpd_id, sol):
        cpd_acc = {}
        metabolite = model.metabolites.get_by_id(cpd_id)
        for r in metabolite.reactions:
            rxn_s = r.metabolites[metabolite]
            v = sol.fluxes[r.id]
            cmps = frozenset(r.compartments)
            if v != 0:
                if cmps not in cpd_acc:
                    cpd_acc[cmps] = 0
                cpd_acc[cmps] += v * rxn_s
        return cpd_acc

    @staticmethod
    def _cpd_acc_to_a_r_t(cpd_acc):
        """
        parse cpd acc to acido, rhoda, total
        """
        return cpd_acc.get(frozenset({'e0'}), 0), cpd_acc.get(frozenset({'e0', 'cA'}), 0), cpd_acc.get(frozenset({'e0', 'cR'}), 0)
