import pandas
import yaml
import subprocess
from datetime import date



input = snakemake.input
output = snakemake.output
params = snakemake.params

home_input_dict = {name:file for name,file in zip(params.home_names,
                                                  input.home)}

with open(input.report) as report_yaml:
    report_dict = yaml.load(report_yaml)

home_pair, home = params.home_list

home_intro_list = []
home_file_list = []
home_stat_list = []


if params.bowtie_options != '':
    bowtie = ''.join(('The following bowtie options were used:\n    ',
                      params.bowtie_options, '\n\n'))
else:
    bowtie = ''

if home_pair in params.pairing_df['id'].values:
    home_df = params.pairing_df[params.pairing_df['id']==home_pair]
    sample_list = home_df.ix[:,'forward':'reverse'].values[0]

    home_id_dict = {'ID':home_pair, 'forward':sample_list[0],
                    'reverse':sample_list[1]}

    intro_str = report_dict['intro_home_tagmap'].format(**home_id_dict)
    home_intro_list.append(intro_str)
    home_dict = {'insert_home':home_input_dict['insert']}
    if sample_list[0] != '-':
        home_dict.update({'pars_structure_fwd':home_input_dict['structure_fwd'],
                          'pars_stats_fwd':home_input_dict['parsed_fwd'],
                          'map_stats_fwd':home_input_dict['maplog_fwd'],
                          'map_fwd':home_input_dict['map_fwd'],
                          'markdup_stats_fwd':home_input_dict['dup_fwd'],
                          'sorted_fwd':home_input_dict['sort_fwd']})
    else:
        home_dict.update({'pars_structure_fwd':"",
                          'pars_stats_fwd':"",
                          'map_stats_fwd':"",
                          'map_fwd':"",
                          'markdup_stats_fwd':"",
                          'sorted_fwd':""})
    if sample_list[1] != '-':
        home_dict.update({'pars_structure_rev':home_input_dict['structure_rev'],
                          'pars_stats_rev':home_input_dict['parsed_rev'],
                          'map_stats_rev':home_input_dict['maplog_rev'],
                          'map_rev':home_input_dict['map_rev'],
                          'markdup_stats_rev':home_input_dict['dup_rev'],
                          'sorted_rev':home_input_dict['sort_rev']})
    else:
        home_dict.update({'pars_structure_rev':"",
                          'pars_stats_rev':"",
                          'map_stats_rev':"",
                          'map_rev':"",
                          'markdup_stats_rev':"",
                          'sorted_rev':""})
    file_str = report_dict['files_home_tagmap'].format(**home_dict)
    home_file_list.append(file_str)

    stats_dict = report_dict['stats_home_tagmap']
    body_list = [body.format(**home_id_dict, **home_dict,
                             bowtie_options=bowtie)
                 if key.startswith('pattern')
                 else body for key,body in stats_dict.items()]
    stat_str = '\n'.join(body_list)
    home_stat_list.append(stat_str)
elif 'insert_arm' in home_input_dict:
    intro_str = report_dict['stats_arms']
    home_intro_list.append(report_dict['intro_arms'])

    file_fmt = report_dict['files_arms']
    file_str = file_fmt.format(arm_region = home_input_dict['insert_arm'])
    home_file_list.append(file_str)
    home_stat_list.append(report_dict['stats_arms'])

this_df = params.pairing_df[params.pairing_df['id']==params.pair]
sample_list = this_df.ix[:,'forward':'reverse'].values[0]



form_dict = {'date': date.today(),
             'ID': params.pair,
             'home_ID': home_pair,
             'forward': sample_list[0],
             'reverse': sample_list[1],
             'home_type': params.home_type,
             'home_intro': '\n'.join(home_intro_list),
             'pars_structure_fwd': input.structure[0],
             'pars_structure_rev': input.structure[1],
             'pars_stats_fwd': input.parsed[0],
             'pars_stats_rev': input.parsed[1],
             'map_stats_fwd': input.maplog[0],
             'map_stats_rev': input.maplog[1],
             'map_fwd': input.map[0],
             'map_rev': input.map[1],
             'markdup_stats_fwd': input.dup[0],
             'markdup_stats_rev': input.dup[1],
             'sorted_fwd': input.sort[0],
             'sorted_rev': input.sort[1],
             'home_exp_files': '\n'.join(home_file_list),
             'random_dist': input.random,
             'mapped_read_dist_fwd': input.mapped_read[0],
             'mapped_read_dist_rev': input.mapped_read[1],
             'sorted_read_dist_fwd': input.sorted_read[0],
             'sorted_read_dist_rev': input.sorted_read[1],
             'insert_dist': input.insert_dist,
             'insert': input.insert,
             'bowtie_options': bowtie,
             'genome_fai': input.fai,
             'home_stats': '\n'.join(home_stat_list)}

body_list = [body.format(**form_dict) if key.startswith('pattern')
             else body for key,body in report_dict['body'].items()]
with open(output.rmd, 'w') as f_out:
    print('\n'.join(body_list), file=f_out)

subprocess.run([input.make_r, output.rmd, output.html])
