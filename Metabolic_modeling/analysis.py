from tqdm import tqdm
import pandas as pd
import cobra
from modelseedpy.core.msatpcorrection import MSATPCorrection
from modelseedpy.core.msgenome import normalize_role
from cobra.core import Reaction
from cobra.flux_analysis import pfba

def calc_max_ATPM(model, exp_d):
    model.reactions.bio1.lower_bound = 0
    model.reactions.ATPM_c0.lower_bound = 0
    model.reactions.EX_cpd00209_e0.lower_bound = exp_d['no3_uptake']
    model.reactions.EX_cpd00029_e0.lower_bound = exp_d['ac_uptake']

    for k in {(list(model.reactions.get_by_id(x[0]).metabolites)[0].name, x[1]) for x in model.medium.items()}:
        print(k)
    model.reactions.bio1.lower_bound = exp_d['growth_OD'] * exp_d['OD_coeff']
    model.objective = 'ATPM_c0'
    sol_max_atpm = model.optimize()
    exp_d['predicted_ATPM'] = sol_max_atpm.fluxes['ATPM_c0']
    print('MAX ATPM', exp_d['predicted_ATPM'])
    model.reactions.bio1.lower_bound = 0
    model.reactions.ATPM_c0.lower_bound = exp_d['predicted_ATPM']
    model.objective = 'bio1'
    sol_max_wt = model.optimize()
    exp_d['predicted_growht_gDW'] = sol_max_wt.fluxes['bio1']
    print('biomass gDW', sol_max_wt.fluxes['bio1'], 'OD600', sol_max_wt.fluxes['bio1'] / exp_d['OD_coeff'])
    return sol_max_atpm, sol_max_wt

SEED_OTHER = [
    'cpd00010_c0',  # coa
    'cpd00118_c0',  # putrs
    'cpd00264_c0',  # spme
    'cpd00056_c0',  # TPP
    'cpd00016_c0',  # pydx
    'cpd00003_c0',  # nad
    'cpd00015_c0',  # fad
    'cpd00042_c0',  # GSH
    'cpd03736_c0',  # Lipid IV
    'cpd15352_c0',  # DMQ
    'cpd15793_c0',  # clp
    'cpd00345_c0',  # 5-Methyltetrahydrofolate 
    'cpd00087_c0',  # THF
    'cpd00201_c0',  # 10-Formyltetrahydrofolate 
    
    'cpd00356_c0',  # dCTP
    'cpd00357_c0',  # TTP
    'cpd00241_c0',  # dGTP
    'cpd00115_c0',  # dATP
    'cpd00038_c0',  # GTP
    'cpd00052_c0',  # CTP
    'cpd00062_c0',  # UTP
    'cpd00002_c0',  # ATP
]

SEED_AA_FLAVOR = {
            'pyr': 'cpd00020_c0',
            'leu': 'cpd00107_c0',
            'lys': 'cpd00039_c0',
            'his': 'cpd00119_c0',
            'ile': 'cpd00322_c0',
            'thr': 'cpd00161_c0',
            'trp': 'cpd00065_c0',
            'tyr': 'cpd00069_c0',
            'ser': 'cpd00054_c0',
            'met': 'cpd00060_c0',
            'cys': 'cpd00084_c0',
            'arg': 'cpd00051_c0',
            'asn': 'cpd00132_c0',
            'asp': 'cpd00041_c0',
            'ala': 'cpd00035_c0',
            'gln': 'cpd00053_c0',
            'glu': 'cpd00023_c0',
            'gly': 'cpd00033_c0',
            'val': 'cpd00156_c0',
            'pro': 'cpd00129_c0',
            'phe': 'cpd00066_c0',
        }

class ProfilerAA:

    def __init__(self, model):
        self.model = model
        #self.amino_acids = BIGG_AA_FLAVOR
        #self.others = BIGG_FLAVOR_2
        self.test_reactions = {}

    def build_synth_metabolite_test(self, metabolite):
        if metabolite.id in self.model.metabolites:
            m = self.model.metabolites.get_by_id(metabolite.id)
            rxn_test = Reaction(f'test_{m.id}', f'Test {m.id} [{m.name}]', 'TEST', 0, 0)
            rxn_test.add_metabolites({
                m: -1
            })
            return rxn_test
        else:
            return None

    def build_tests(self):
        for aa, cpd_id in self.amino_acids.items():
            cpd = self.model.metabolites.get_by_id(cpd_id)
            rxn_test = Reaction(f'test_{aa}', f'Test {cpd.id} [{cpd.name}]', 'TEST', 0, 0)
            rxn_test.add_metabolites({
                cpd: -1
            })

            if rxn_test.id not in self.model.reactions:
                self.test_reactions[aa] = rxn_test
        self.model.add_reactions(list(self.test_reactions.values()))

    def profile_genome(self, genome_id):
        result = {}
        for test_id, r in self.test_reactions.items():
            self.model.objective = r.id
            model_reaction = self.reactions.get_by_id(r.id)
            model_reaction.lower_bound = 0
            model_reaction.upper_bound = 1
            solution = pfba(self.model)
            obj = solution.fluxes[r.id]
            if obj > 0 and solution.status == 'optimal':
                result[test_id] = solution
            else:
                result[r.id] = None
            r.upper_bound = 0

        return result

def _isolate_test(model_orig):
    model = cobra.io.from_json(cobra.io.to_json(model_orig))
    model.reactions.EX_cpd00209_e0.lower_bound = -1000
    model.reactions.EX_cpd00029_e0.lower_bound = -1000
    model.reactions.ATPM_c0.lower_bound = 0
    prof = ProfilerAA(model)

    tests = []
    for v in SEED_AA_FLAVOR.values():
        m = model.metabolites.get_by_id(v)
        rxn = prof.build_synth_metabolite_test(m)
        if rxn:
            tests.append(rxn)
    for v in SEED_OTHER:
        m = model.metabolites.get_by_id(v)
        rxn = prof.build_synth_metabolite_test(m)
        if rxn:
            tests.append(rxn)
    model.add_reactions(tests)

    solutions = {}
    for t in tqdm(tests):
        rxn = model.reactions.get_by_id(t.id)
        rxn.upper_bound = 1
        model.objective = rxn.id
        sol = cobra.flux_analysis.pfba(model)
        solutions[rxn.id] = sol
        rxn.upper_bound = 0

    return solutions, tests

def _report(solutions_acido, solutions_rhoda, tests):
    rows = []
    for r in tests:
        if r.id.startswith('test_'):
            
            m = list(r.metabolites)[0]
            sol_acido = solutions_acido[r.id]
            sol_rhoda = solutions_rhoda[r.id]
            atp_acido = sol_acido.fluxes['rxn08173_c0']
            atp_rhoda = sol_rhoda.fluxes['rxn08173_c0']
            
            ac_acido = sol_acido.fluxes['EX_cpd00029_e0']
            ac_rhoda = sol_rhoda.fluxes['EX_cpd00029_e0']
            nitrate_acido = sol_acido.fluxes['EX_cpd00209_e0']
            nitrate_rhoda = sol_rhoda.fluxes['EX_cpd00209_e0']
            
            data = [m.id, m.name, sol_acido.fluxes[r.id], sol_rhoda.fluxes[r.id], 
                    atp_acido, atp_rhoda,
                    ac_acido, ac_rhoda,
                    nitrate_acido, nitrate_rhoda,
                   ]
            #print(data)
            data = {
                'cpd_id': m.id,
                'cpd_name': m.name,
                '3H11': sol_acido.fluxes[r.id],
                'R12': sol_rhoda.fluxes[r.id],
                '3H11 (ATP)': atp_acido,
                'R12 (ATP)': atp_rhoda,
                '3H11 (acetate)': ac_acido,
                'R12 (acetate)': ac_rhoda,
                '3H11 (nitrate)': nitrate_acido,
                'R12 (nitrate)': nitrate_rhoda,
            }
            rows.append(data)
    df = pd.DataFrame(rows)
    df = df.set_index('cpd_id')
    return df

def model_comm_io2(model_comm, solution_comm_wt):
    for e in model_comm.exchanges:
        v = solution_comm_wt.fluxes[e.id]
        if v != 0:
            m = list(e.metabolites)[0]
            print(m.id, m.name, v)
            total = {}
            for rxn in m.reactions:
                rxn_s = rxn.metabolites[m]
                rxn_v = solution_comm_wt.fluxes[rxn.id]
                if rxn_v != 0:
                    cmps = frozenset(rxn.compartments)
                    if cmps not in total:
                        total[cmps] = 0
                    total[cmps] += rxn_v * rxn_s
            total_v = 0
            for cmps in total:
                cmp_str = ';'.join(cmps)
                
                if not 'e0' == cmp_str:
                    total_v += total[cmps]
                    print('\t', cmp_str, total[cmps])
            print('\t', f'total {m.name}', total_v)
            
def model_comm_io(model_comm, solution_comm_wt, exclude=None):
    if exclude is None:
        exclude = set()
    for m in model_comm.metabolites:
        if m.id not in exclude and m.compartment == 'e0':
            total = {}
            for rxn in m.reactions:
                rxn_s = rxn.metabolites[m]
                rxn_v = solution_comm_wt.fluxes[rxn.id]
                if rxn_v != 0:
                    cmps = frozenset(rxn.compartments)
                    if cmps not in total:
                        total[cmps] = 0
                    total[cmps] += rxn_v * rxn_s
            total_v = 0
            print_total = False
            for cmps in total:
                cmp_str = ';'.join(cmps)
                if not 'e0' == cmp_str:
                    total_v += total[cmps]
                    if round(total[cmps], 9) != 0:
                        print_total = True
                        print(m.id, m.name, cmp_str, total[cmps])
            if print_total:
                print('\t', f'total {m.name}', round(total_v, 10))