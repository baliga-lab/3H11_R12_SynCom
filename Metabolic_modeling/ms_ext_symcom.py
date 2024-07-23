from modelseedpy import MSBuilder, MSATPCorrection, MSMedia, MSGapfill


MEDIA_GENOME_SCALE = medium_gapfill = {
    'EX_cpd00029_e0': 20.0,
 'EX_cpd00013_e0': 100.0,
 'EX_cpd00001_e0': 100.0,
 'EX_cpd00218_e0': 0.0, #niacin
 'EX_cpd00220_e0': 100.0,
 'EX_cpd00305_e0': 100.0,
 'EX_cpd00393_e0': 100.0,
 'EX_cpd03424_e0': 100.0,
 'EX_cpd00443_e0': 0.0, #ABEE
 'EX_cpd00644_e0': 0.0002281,
 'EX_cpd00263_e0': 100.0,
 'EX_cpd00048_e0': 100.0,
 'EX_cpd00009_e0': 100.0,
 'EX_cpd00242_e0': 29.759425,
 'EX_cpd00205_e0': 1.3415688,
 'EX_cpd00063_e0': 100.0,
 'EX_cpd00971_e0': 34.9324073,
 'EX_cpd00099_e0': 100.0,
 'EX_cpd00254_e0': 100.0,
 'EX_cpd00030_e0': 100.0,
 'EX_cpd00058_e0': 100.0,
 'EX_cpd00034_e0': 100.0,
 'EX_cpd10515_e0': 100.0,
 'EX_cpd00149_e0': 100.0,
 'EX_cpd00244_e0': 100.0,
 'EX_cpd11574_e0': 100.0,
 'EX_cpd15574_e0': 100.0,
 'EX_cpd00067_e0': 100.0,
 'EX_cpd00209_e0': 10.0}


def load_genome_from_annotation_file(filename):
    from pandas import read_csv
    from modelseedpy.core.msgenome import MSGenome, MSFeature

    df = read_csv(filename, sep='\t', index_col=0)
    features = {}
    for feature_id, d in df.iterrows():
        gene_id = d['gene_id']
        feature = MSFeature(gene_id, '')
        for rast_feature in d['RAST'].split('; '):
            feature.add_ontology_term('RAST', rast_feature)
        if feature.id not in features:
            features[feature.id] = feature
        else:
            raise ValueError('duplicate gene ID')
    genome = MSGenome()
    genome.add_features(list(features.values()))

    return genome


from cobra.core import Metabolite, Reaction
from modelseedpy import MSBuilder

def _add_leucine_import_abc(model):
    if 'cpd00107_e0' not in model.metabolites:
        model.add_metabolites([Metabolite('cpd00107_e0', 'C6H13NO2', 'L-Leucine [e0]', 0, 'e0')])
        print(f'[{model.id}] added L-Leucine [e0]')
    if 'rxn05161_c0' not in model.reactions:
        bcaa_leu_abc_transporter = Reaction('rxn05161_c0', 'ABC Transporter (L-Leucine) [c0]', 'Transport', 0, 1000)
        bcaa_leu_abc_transporter.add_metabolites({
            model.metabolites.cpd00001_c0: -1,
            model.metabolites.cpd00002_c0: -1,
            model.metabolites.cpd00107_e0: -1,
            model.metabolites.cpd00107_c0: 1,
            model.metabolites.cpd00008_c0: 1,
            model.metabolites.cpd00009_c0: 1,
            model.metabolites.cpd00067_c0: 1,    
        })
        model.add_reactions([bcaa_leu_abc_transporter])
        print(f'[{model.id}] added ABC Transporter (L-Leucine) [c0]')
    MSBuilder.add_exchanges_to_model(model)
def _add_leucine_export(model):
    if 'cpd00107_e0' not in model.metabolites:
        model.add_metabolites([Metabolite('cpd00107_e0', 'C6H13NO2', 'L-Leucine [e0]', 0, 'e0')])
        print(f'[{model.id}] added L-Leucine [e0]')
    if 'LeuE_c0' not in model.reactions:
        leuE_transporter = Reaction('LeuE_c0', 'LeuE [c0]', 'Transport', -1000, 1000)
        leuE_transporter.add_metabolites({
            model.metabolites.cpd00107_e0: -1,
            model.metabolites.cpd00107_c0: 1, 
            model.metabolites.cpd00067_e0: 1,
            model.metabolites.cpd00067_c0: -1, 
        })
        model.add_reactions([leuE_transporter])
        print(f'[{model.id}] added LeuE')
    MSBuilder.add_exchanges_to_model(model)


def _setup_fix(model):
    delete = [
        'rxn09003_c0', 'rxn09001_c0', 'rxn05627_c0',
        'rxn08062_c0', # Acetate [e0] + Na+ [e0] --> Acetate [c0] + Na+ [c0]
        'rxn05298_c0', # L-Glutamate [e0] + Na+ [e0] <=> L-Glutamate [c0] + Na+ [c0]
        'rxn05313_c0', # Phosphate [e0] + 3 Na+ [e0] --> Phosphate [c0] + 3 Na+ [c0]
    ]
    _to_delete = []
    for i in delete:
        if i in model.reactions:
            _to_delete.append(i)
    model.reactions -= _to_delete
    
    # H+ [e0] + Nitrate [e0] --> H+ [c0] + Nitrate [c0]
    if 'rxn05627_c0' in model.reactions:
        model.reactions.rxn05627_c0.upper_bound = 0
        
    # 2 Cytochrome c3+ [c0] + D-Lactate [c0] --> Pyruvate [c0] + 2 H+ [c0] + 2 Cytochrome c2+ [c0]
    if 'rxn00146_c0' in model.reactions:
        model.reactions.rxn00146_c0.lower_bound = 0

    # Phosphoenolpyruvate [c0] + GLUM [e0] --> Pyruvate [c0] + D-Glucosamine phosphate [c0]
    if 'rxn05569_c0' in model.reactions:
        model.reactions.rxn05569_c0.lower_bound = 0
    
    # FAD [c0] + H+ [c0] + Isovaleryl-CoA [c0] --> FADH2 [c0] + Dimethylacryloyl-CoA [c0]
    if 'rxn02866_c0' in model.reactions:
        model.reactions.rxn02866_c0.lower_bound = 0
        model.reactions.rxn02866_c0.upper_bound = 1000


class SynComStudy:

    def __init__(self, template):
        self.rast = None
        self.kbase = None
        self.genome_acido = None
        self.genome_rhoda = None
        self.template = template

    def build_isolates(self):
        model_acido = self.build_isolate('3H11', self.genome_acido)
        model_rhoda = self.build_isolate('R12', self.genome_rhoda)
        _setup_fix(model_acido)
        _setup_fix(model_rhoda)
        _add_leucine_import_abc(model_rhoda)
        _add_leucine_export(model_acido)
        model_rhoda.reactions.EX_cpd00107_e0.lower_bound = 0

        model_rhoda.reactions.dnr00001_c0.add_metabolites({
            model_rhoda.metabolites.cpd00067_c0: 2,
            model_rhoda.metabolites.cpd00067_e0: -2
        })
        print('R12', model_rhoda.reactions.dnr00001_c0.build_reaction_string(True))
        print('3H11', model_acido.reactions.dnr00001_c0.build_reaction_string(True))

        return model_acido, model_rhoda

    def build_isolate(self, model_id, genome):
        model = self.build_model(model_id, genome)
        model_gapfill = model.copy()

        atp_tests = self.atp_correction(model)
        gap_fill = self.gap_fill(model, atp_tests)

        rxn_new, ex_new = self._integrate_solution(model_gapfill, gap_fill)
        print(f'gapfill reactions {len(rxn_new)} new exchanges {len(ex_new)}')
        model_gapfill.medium = {k: v for k, v in MEDIA_GENOME_SCALE.items() if k in model_gapfill.reactions}

        return model_gapfill

    @staticmethod
    def _add_atpm(model):
        from cobra.core import Reaction
        if 'ATPM_c0' not in model.reactions:
            atpm = Reaction(f'ATPM_c0', f'ATPM', 'ATPM', 0, 1000)
            atpm.add_metabolites({
                model.metabolites.cpd00001_c0: -1,
                model.metabolites.cpd00002_c0: -1,
                model.metabolites.cpd00008_c0: 1,
                model.metabolites.cpd00009_c0: 1,
                model.metabolites.cpd00067_c0: 1,
            })
            model.add_reactions([atpm])

    def build_model(self, model_id, genome):
        """
        build base model
        :param model_id:
        :param genome:
        :return:
        """
        b = MSBuilder(genome, self.template, model_id)
        model_base = b.build(model_id, annotate_with_rast=False)
        SynComStudy._add_atpm(model_base)
        return model_base

    def atp_correction(self, model):
        media = MSMedia.from_dict({
            'cpd00001': 1000,
            'cpd00067': 1000,
            'cpd00209': 1,
            'cpd00029': 1,
        })
        media.id = 'nitrate'
        media.name = 'nitrate'
        media.get_media_constraints()
        atp_correction = MSATPCorrection(model, self.template, [media],
                                         'c0', atp_hydrolysis_id='ATPM_c0',
                                         load_default_medias=False)
        tests = atp_correction.run_atp_correction()

        return tests

    def _get_solution(self, gf, sol, model):
        res = {}
        for rxn in gf.gfmodel.reactions:
            v = round(sol.fluxes[rxn.id], 9)
            if v != 0:
                if rxn.id not in model.reactions:
                    if rxn.id[:-1] in self.template.reactions:
                        # print(v, rxn.id, rxn.build_reaction_string(True))
                        lb = -1000
                        ub = 1000
                        if v > 0:
                            lb = 0
                        elif v < 0:
                            ub = 0
                        res[rxn.id[:-1]] = (lb, ub)
        return res

    def _integrate_solution(self, model, gap_fill_solution):
        added_reactions = []
        for rxn_id, (lb, ub) in gap_fill_solution.items():
            template_reaction = self.template.reactions.get_by_id(rxn_id)
            model_reaction = template_reaction.to_reaction(model)
            model_reaction.lower_bound = lb
            model_reaction.upper_bound = ub
            added_reactions.append(model_reaction)
        model.add_reactions(added_reactions)
        add_exchanges = MSBuilder.add_exchanges_to_model(model)

        return added_reactions, add_exchanges

    def gap_fill(self, model, tests, min_biomass=0.02):
        model.objective = 'bio1'
        gapfill = MSGapfill(model,
                            default_gapfill_templates=[self.template],
                            test_conditions=tests,
                            default_target='bio1')

        gapfill.gfmodel.reactions.bio1.lower_bound = min_biomass
        gapfill.gfmodel.medium = MEDIA_GENOME_SCALE

        gapfill_fba = gapfill.gfmodel.optimize()

        gapfill_solution = self._get_solution(gapfill, gapfill_fba, model)

        return gapfill_solution


def build_genome_annotation_file(genome, gene_id_remap, starts_with_filter='FW510-R12',
                                 filename_psort='psortb_r12_results.tsv',
                                 filename_out='annotation_rhoda.tsv'):
    """
    :param genome: rhoda or acido genome
    :param gene_id_remap: dict mapping feature_id to id from expression data
    :param starts_with_filter: FW510-R12 for rhoda GW101-3H11 for acido
    :param filename_psort: psort output file
    :param filename_out: output annotation_file
    :return:
    """
    psort_data = {}
    with open(filename_psort, 'r') as fh:
        h = fh.readline()
        l = fh.readline()
        header = {v: i for i, v in enumerate(h.split('\t'))}
        while l:
            _p = [x.strip() for x in l.split('\t')]
            psort_data[_p[header['SeqID']]] = [
                _p[header['Final_Localization']],
                _p[header['Final_Score']],
                _p[header['Final_Localization_Details']],
                _p[header['Secondary_Localization']],
            ]
            l = fh.readline()
    remap = {}
    for k, v in gene_id_remap.items():
        if k.startswith(starts_with_filter):
            if v not in remap:
                remap[v] = k
            else:
                print(k)
    with open(filename_out, 'w') as fh:
        d = [
            'feature_id', 'gene_id', 'RAST',
            'psort_loc', 'psort_loc_score', 'psort_loc_details', 'psort_loc_sec']
        fh.write('\t'.join(d) + '\n')
        for f in genome.features:
            d = [
                f.id,
                remap[f.id],
                '; '.join(f.ontology_terms.get('RAST', [])),
                psort_data[f.id][0],
                psort_data[f.id][1],
                psort_data[f.id][2],
                psort_data[f.id][3]
            ]
            fh.write('\t'.join(d) + '\n')


class GSPComBuilder:
    """
    Assembles Community model for denitrification
    """

    def __init__(self, model_rhoda, model_acido, media_com):
        """

        :param model_rhoda: Model
        :param model_acido: Model
        :param media_com: Media
        """
        self.models = {model_rhoda: "R", model_acido: "A"}
        self.media_com = media_com

    def a(self, model_com):
        ## Replace to cytochromes
        model_com.reactions.rxn01806_cA.add_metabolites(
            {
                model_com.metabolites.cpd00004_cA: +1,
                model_com.metabolites.cpd00003_cA: -1,
                model_com.metabolites.cpd18072_cA: +1,
                model_com.metabolites.cpd18074_cA: -1,
            }
        )
        model_com.reactions.rxn01806_cR.add_metabolites(
            {
                model_com.metabolites.cpd00004_cR: +1,
                model_com.metabolites.cpd00003_cR: -1,
                model_com.metabolites.cpd18072_cR: +1,
                model_com.metabolites.cpd18074_cR: -1,
            }
        )

    def external_denitri(self, model):
        from cobra.core import Metabolite, Reaction

        print("external denitri")
        delete = {
            "rxn10577_cR",
            "rxn10577_cA",
            "rxn08966_cA",
            "rxn08966_eR",
            "rxn09008_cR",
            "rxn09008_cA",
            "rxn05890_cR",
            "rxn11932_cA",
            "rxn11932_cR",
            "rxn05627_cR",
            "rxn05627_cA",
            "rxn05625_cA",
        }
        model.remove_reactions(delete)
        model.add_metabolites([Metabolite("cpd00418_e0", "NO", "NO_e0")])
        rxnR = Reaction("rxn09004_cR", "rxn09004_cR", "", -1000, 1000)
        rxnR.add_metabolites(
            {
                model.metabolites.cpd00209_cR: 1,
                model.metabolites.cpd00075_cR: -1,
                model.metabolites.cpd00209_e0: -1,
                model.metabolites.cpd00075_e0: 1,
            }
        )
        rxnA = Reaction("rxn09004_cA", "rxn09004_cA", "", -1000, 1000)
        rxnA.add_metabolites(
            {
                model.metabolites.cpd00209_cA: 1,
                model.metabolites.cpd00075_cA: -1,
                model.metabolites.cpd00209_e0: -1,
                model.metabolites.cpd00075_e0: 1,
            }
        )
        ex = Reaction("EX_cpd00418_e0", "EX_cpd00418_e0", "", 0, 1000)
        ex.add_metabolites({model.metabolites.cpd00418_e0: -1})
        model.add_reactions([rxnR, rxnA, ex])
        move = {
            "cpd00528_cR": ["cpd00528_e0", None],
            "cpd00528_cA": ["cpd00528_e0", None],
            "cpd00659_cR": ["cpd00659_e0", None],
            "cpd00659_cA": ["cpd00659_e0", None],
            "cpd00418_cA": ["cpd00418_e0", None],
            "cpd00418_cR": ["cpd00418_e0", None],
            "cpd00075_cR": ["cpd00075_e0", {"rxn14428_cR"}],
            # 'cpd00075_cR': ['cpd00075_e0', {'rxn09001_cR', 'rxn14427_cR', 'rxn09003_cR', 'rxn14428_cR'}],
        }

        for src_cpd_id in move:
            dst_cpd_id = move[src_cpd_id][0]
            if dst_cpd_id not in model.metabolites:
                print("!")
                break
            dst_cpd = model.metabolites.get_by_id(dst_cpd_id)
            src_cpd = model.metabolites.get_by_id(src_cpd_id)
            target_rxn_ids = move[src_cpd_id][1]
            for rxn in src_cpd.reactions:
                if target_rxn_ids is None or rxn.id in target_rxn_ids:
                    rxn.add_metabolites(
                        {
                            src_cpd: -1 * rxn.metabolites[src_cpd],
                            dst_cpd: rxn.metabolites[src_cpd],
                        }
                    )
        return model

    def b(self, model_com):
        ## Replace to cytochromes
        model_com.reactions.rxn01806_cA.add_metabolites(
            {
                model_com.metabolites.cpd00004_cA: +1,
                model_com.metabolites.cpd00003_cA: -1,
                model_com.metabolites.cpd15560_cA: +1,
                model_com.metabolites.cpd15561_cA: -1,
            }
        )
        model_com.reactions.rxn01806_cR.add_metabolites(
            {
                model_com.metabolites.cpd00004_cR: +1,
                model_com.metabolites.cpd00003_cR: -1,
                model_com.metabolites.cpd15560_cR: +1,
                model_com.metabolites.cpd15561_cR: -1,
            }
        )

    def build(self):
        from modelseedpy import MSBuilder

        ex_ids = set(self.media_com.get_media_constraints().keys())
        ex_ids.add("cpd00075_e0")
        ex = {}
        m_pointer = {}

        from cobra.core import Model, Reaction, Metabolite

        model_com = Model("com")
        for model, token in self.models.items():
            metabolites = []
            for m in model.metabolites:
                if m.id in ex_ids:
                    # print(m)
                    if m.id not in ex:
                        m_copy = Metabolite(
                            m.id, m.formula, m.name, m.charge, m.compartment
                        )
                        m_pointer[m.id] = m.id
                        metabolites.append(m_copy)
                else:
                    m_id_copy = m.id[:-1] + token
                    # print(m_id_copy, m.id)
                    m_copy = Metabolite(
                        m_id_copy,
                        m.formula,
                        m.name[:-1] + token,
                        m.charge,
                        m.compartment[:-1] + token,
                    )
                    m_pointer[m.id] = m_copy.id
                    metabolites.append(m_copy)
            model_com.add_metabolites(metabolites)
            reactions = []
            for r in model.reactions:
                if not r.id.startswith("EX_"):
                    if r.id.startswith("bio"):
                        r_copy_id = r.id + "_" + token
                    else:
                        r_copy_id = r.id[:-1] + token
                    r_copy = Reaction(
                        r_copy_id,
                        r.name[:-1] + token,
                        r.subsystem,
                        r.lower_bound,
                        r.upper_bound,
                    )
                    r_copy_m = dict(
                        map(
                            lambda x: (
                                model_com.metabolites.get_by_id(m_pointer[x[0].id]),
                                x[1],
                            ),
                            r.metabolites.items(),
                        )
                    )
                    r_copy.add_metabolites(r_copy_m)
                    reactions.append(r_copy)
                    # print(r_copy.name)
                else:
                    # print(r)
                    pass
            model_com.add_reactions(reactions)
            MSBuilder.add_exchanges_to_model(model_com, 'e' + token)

        MSBuilder.add_exchanges_to_model(model_com, 'e0')
        print('Exchanges ', len(model_com.exchanges), 'e0')

        r_bio_sum = Reaction("bio1", "bio_com", "", 0, 1000)
        r_bio_sum.add_metabolites(
            {
                model_com.metabolites.cpd11416_cA: -0.4,
                model_com.metabolites.cpd11416_cR: -0.6,
            }
        )
        rxn09008_cR = Reaction("rxn09008_cR", "rxn09008_cR", "", -1000, 1000)
        rxn09008_cR.add_metabolites(
            {
                model_com.metabolites.cpd00418_cR: -1,
                model_com.metabolites.cpd00418_eA: 1,
            }
        )
        model_com.add_reactions([r_bio_sum, rxn09008_cR])
        model_com.objective = "bio1"
        model_com.reactions.bio1.annotation = {"sbo": "SBO:0000629"}
        model_com.reactions.bio1_A.annotation = {"sbo": "SBO:0000629"}
        model_com.reactions.bio1_R.annotation = {"sbo": "SBO:0000629"}
        for r in model_com.exchanges:
            r.lower_bound = 0
            r.upper_bound = 1000
        for k, (lb, ub) in self.media_com.get_media_constraints().items():
            if "EX_" + k in model_com.reactions:
                r = model_com.reactions.get_by_id("EX_" + k)
                r.lower_bound = lb
                r.upper_bound = ub
        model_com.reactions.EX_cpd00528_e0.upper_bound = 100
        model_com.reactions.EX_cpd00659_e0.upper_bound = 100
        model_com.reactions.rxn05488_cA.add_metabolites(
            {
                model_com.metabolites.cpd00067_e0: +1,
                model_com.metabolites.cpd00067_cA: -1,
            }
        )
        # rxn = Reaction('rxn05625_cR', 'rxn05625_cR', '', -1000, 1000)
        # rxn.add_metabolites({model_com.metabolites.cpd00075_cR: -1, model_com.metabolites.cpd00075_e0: 1})
        # model_com.add_reactions([rxn])

        # model_com.reactions.rxn11937_cR.lower_bound = 0 #-1000 # block nosZ
        return model_com

    def media_fix(self, model_com, media_com):
        mediacompounds = []
        for ex_id, v in model_com.medium.items():
            compound_id = ex_id[:-3][3:]
            print(
                list(
                    filter(
                        lambda x: x["id"] == compound_id + "_e0",
                        model_o["modelcompounds"],
                    )
                )[0]["compound_ref"]
            )
            mediacompounds.append(
                {
                    "compound_ref": "489/6/15/compounds/id/" + compound_id,
                    "concentration": 20,
                    "id": compound_id,
                    "inchikey": "",
                    "maxFlux": v,
                    "minFlux": -100,
                    "name": compound_id,
                    "smiles": "",
                }
            )
        data = media_com.get_data()
        data["mediacompounds"] = mediacompounds
        return data

    def fix_kbase_object_data(self, model_acido, model_o, model1, model2):
        """
        model1: model_com
        model2: kbase_model
        model_o: kbase_model (object)
        """
        model_o["genome_ref"] = model_acido.genome_ref
        model_o["template_ref"] = model_acido.template_ref

        for r in model1.reactions:
            if r.id not in model2.reactions:
                cpd_id = list(r.metabolites)[0]
                if len(r.metabolites) > 1:
                    print(cpd_id.id, r.metabolites)
                else:
                    modelReactionReagents = list(
                        map(
                            lambda x: {
                                "coefficient": x[1],
                                "modelcompound_ref": "~/modelcompounds/id/" + x[0].id,
                            },
                            r.metabolites.items(),
                        )
                    )
                    print(cpd_id.id, modelReactionReagents)

                    compound_id = cpd_id.id
                    mr = {
                        "aliases": [],
                        "dblinks": {},
                        "direction": ">",
                        "edits": {},
                        "gapfill_data": {},
                        "id": "OUT_" + compound_id + "e",
                        "maxforflux": 1000,
                        "maxrevflux": 0,
                        "modelReactionProteins": [],
                        "modelReactionReagents": modelReactionReagents,
                        "modelcompartment_ref": "~/modelcompartments/id/c0",
                        "name": "SINK OF " + compound_id,
                        "numerical_attributes": {},
                        "probability": 0,
                        "protons": 0,
                        "reaction_ref": "~/template/reactions/id/rxn00000_c",
                        "string_attributes": {},
                    }
                    model_o["modelreactions"].append(mr)


def comm(models):
    from cobra.core import Model, Reaction, Metabolite
    model_com = Model('com')
    ex_ids = set()
    for model in models.values():
        ex_ids |= {list(r.metabolites)[0].id for r in model.exchanges}

    ex = {}
    m_pointer = {}
    for token, model in models.items():
        m_mapper = {}
        for m in model.metabolites:
            if m.id in ex_ids:
                # print(m)
                if m.id not in ex:
                    m_copy = Metabolite(m.id, m.formula, m.name, m.charge, m.compartment)
                    m_mapper[m.id] = m_copy
            else:
                m_id_copy = m.id[:-1] + token
                # print(m_id_copy, m.id)
                m_copy = Metabolite(m_id_copy, m.formula, m.name[:-1] + token, m.charge, m.compartment[:-1] + token)
                m_mapper[m.id] = m_copy
        m_pointer[token] = m_mapper

    for token in models:
        model_com.add_metabolites(list(m_pointer[token].values()))

    r_pointer = {}
    for token, model in models.items():
        r_mapper = {}
        for r in model.reactions:
            if not r.id.startswith('EX_'):
                if r.id.startswith('bio'):
                    r_copy_id = r.id + '_' + token
                else:
                    r_copy_id = r.id[:-1] + token
                r_copy = Reaction(r_copy_id, r.name[:-1] + token, r.subsystem, r.lower_bound, r.upper_bound)
                r_copy_m = dict(map(lambda x: (model_com.metabolites.get_by_id(m_pointer[token][x[0].id].id), x[1]),
                                    r.metabolites.items()))
                r_copy.add_metabolites(r_copy_m)
                r_mapper[r_copy.id] = r_copy
            else:
                # print(r)
                pass
        r_pointer[token] = r_mapper

    for token in models:
        model_com.add_reactions(list(r_pointer[token].values()))

    r_bio_sum = Reaction('bio1', 'bio_com', '', 0, 1000)
    r_bio_sum.add_metabolites({
        model_com.metabolites.cpd11416_cA: -0.4,
        model_com.metabolites.cpd11416_cR: -0.6
    })
    model_com.add_reactions([r_bio_sum])

    model_com.objective = 'bio1'
    model_com.reactions.bio1.annotation = {'sbo': 'SBO:0000629'}
    model_com.reactions.bio1_A.annotation = {'sbo': 'SBO:0000629'}
    model_com.reactions.bio1_R.annotation = {'sbo': 'SBO:0000629'}

    from modelseedpy import MSBuilder

    exchanges = MSBuilder.add_exchanges_to_model(model_com)

    return model_com

from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Tuple, Union
from cobra.core import Model, Reaction, Metabolite


class Comm(Model):

    def __init__(self, id_or_model: Union[str, "Model", None] = None, name: Optional[str] = None):
        super().__init__(id_or_model, name)


class CommFactory:

    def __init__(self):
        self.models = {}

    def with_model(self, model, abundance, index):
        if index in self.models:
            raise ValueError(f'Invalid index {index}. Already taken')
        self.models[index] = (model, abundance)
        return self

    def build_comm_biomass(self):
        r_bio_sum = Reaction('bio1', 'bio_com', '', 0, 1000)
        r_bio_sum.add_metabolites({
            model_com.metabolites.cpd11416_cA: -0.4,
            model_com.metabolites.cpd11416_cR: -0.6
        })
        model_com.add_reactions([r_bio_sum])

        model_com.objective = 'bio1'
        model_com.reactions.bio1.annotation = {'sbo': 'SBO:0000629'}
        model_com.reactions.bio1_A.annotation = {'sbo': 'SBO:0000629'}
        model_com.reactions.bio1_R.annotation = {'sbo': 'SBO:0000629'}

    def build(self):
        model = Comm('model')

        m_pointer = {}

        ex_ids = set()
        for model in models.values():
            ex_ids |= {list(r.metabolites)[0].id for r in model.exchanges}

        for token, model in models.items():
            m_mapper = {}
            for m in model.metabolites:
                if m.id in ex_ids:
                    # print(m)
                    if m.id not in ex:
                        m_copy = Metabolite(m.id, m.formula, m.name, m.charge, m.compartment)
                        m_mapper[m.id] = m_copy
                else:
                    m_id_copy = m.id[:-1] + token
                    # print(m_id_copy, m.id)
                    m_copy = Metabolite(m_id_copy, m.formula, m.name[:-1] + token, m.charge, m.compartment[:-1] + token)
                    m_mapper[m.id] = m_copy
            m_pointer[token] = m_mapper

        for token in models:
            model.add_metabolites(list(m_pointer[token].values()))

        r_pointer = {}
        for token, model in models.items():
            r_mapper = {}
            for r in model.reactions:
                if not r.id.startswith('EX_'):
                    if r.id.startswith('bio'):
                        r_copy_id = r.id + '_' + token
                    else:
                        r_copy_id = r.id[:-1] + token
                    r_copy = Reaction(r_copy_id, r.name[:-1] + token, r.subsystem, r.lower_bound, r.upper_bound)
                    r_copy_m = dict(map(lambda x: (model_com.metabolites.get_by_id(m_pointer[token][x[0].id].id), x[1]),
                                        r.metabolites.items()))
                    r_copy.add_metabolites(r_copy_m)
                    r_mapper[r_copy.id] = r_copy
                else:
                    # print(r)
                    pass
            r_pointer[token] = r_mapper

        for token in models:
            model_com.add_reactions(list(r_pointer[token].values()))

        return model

def _build(self):
    model_comm = Comm('model')

    ex_ids = set()
    for token in self.models:
        model = self.models[token][0]
        ex_ids |= {list(r.metabolites)[0].id for r in model.exchanges}

    m_pointer = {}
    ex = {}
    for token in self.models:
        model = self.models[token][0]
        m_mapper = {}
        for m in model.metabolites:
            if m.id in ex_ids:
                # print(m)
                if m.id not in ex:
                    m_copy = Metabolite(m.id, m.formula, m.name, m.charge, m.compartment)
                    m_mapper[m.id] = m_copy
            else:
                m_id_copy = m.id[:-1] + token
                # print(m_id_copy, m.id)
                m_copy = Metabolite(m_id_copy, m.formula, f'{m.name[:-4]}[{token}]', m.charge, m.compartment[:-1] + token)
                m_mapper[m.id] = m_copy
        m_pointer[token] = m_mapper

    for token in self.models:
        model_comm.add_metabolites(list(m_pointer[token].values()))

    r_pointer = {}
    for token in self.models:
        model = self.models[token][0]
        r_mapper = {}
        for r in model.reactions:
            if not r.id.startswith('EX_'):
                if r.id.startswith('bio'):
                    r_copy_id = r.id + '_' + token
                else:
                    r_copy_id = r.id[:-1] + token
                r_copy = Reaction(r_copy_id, f'{r.name[:-4]}[{token}]', r.subsystem, r.lower_bound, r.upper_bound)
                # r_copy.gene_reaction_rule = r.gene_reaction_rule
                r_copy_m =  {model_comm.metabolites.get_by_id(m_pointer[token][m.id].id): v for m, v in r.metabolites.items()}
                r_copy.add_metabolites(r_copy_m)
                r_mapper[r_copy.id] = r_copy
            else:
                pass
        r_pointer[token] = r_mapper

    for token in self.models:
        model_comm.add_reactions(list(r_pointer[token].values()))

    MSBuilder.add_exchanges_to_model(model_comm)

    return model_comm