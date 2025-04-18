# this is bad practice but i'm doing it anyway
import warnings
warnings.filterwarnings("ignore")

import cPickle
import train_deepARG
import make_XY
import os
import json

# convert alignments.tsv.json and metadata_LS.pkl to usable alignments dictionary
def process_alignments(alignments, mdfile, ids_to_drop=None):
    try:
        ids_to_drop = set(ids_to_drop)
    except:
        ids_to_drop = set()

    # make dictionary mapping IDs to full features string
    mds = cPickle.load(open(mdfile))
    id2md = dict()
    for md in mds['features']:
        id = md.split('|', 1)[0]
        if id not in ids_to_drop:
            id2md[id] = md

    print('ids:', len(id2md.keys()))
    
    # make new alignments
    new_algns = dict()
    seqs_skipped = 0
    algns_skipped = 0
    for k,v in alignments.iteritems():
        # new key: feature string
        # ensure we can find the ID and the sequence is from uniprot
        try:
            new_k = id2md[k].split('|')
            assert new_k[2] == 'UNIPROT'
            # put label at index 2 for make_xy2; the "column" this removes contains no information
            new_k = '|'.join([new_k[0]] + new_k[2:])
        except:
            seqs_skipped += 1
            continue
        new_v = dict()

        # populate new value; v is a list, each e is a dict
        for e in v:
            try:
            # skip alignments against other uniprot sequences + unknown ARDB/CARD sequences
                assert e['dbsubject'] != 'UNIPROT'
                md = id2md[e['subject']]
            except:
                algns_skipped += 1
                continue
            score = float(e['BitScore'])
            new_v[md] = score
        
        # add to final dictionary
        new_algns[new_k] = new_v

    print('samples found:', len(new_algns.keys()))
    print('seqs skipped:', seqs_skipped)
    print('algns skipped:', algns_skipped)

    return new_algns

# run training and save model
def train(dtpath='../../data', version='new', 
          algnfile = 'scripts/db/alignment.tsv.json',
          mdprefix = 'model/v1/metadata_',
          mode='LS',
          process=True,
          ids_to_drop=None):
    os.system("mkdir -p "+dtpath+"/model/"+version)
    
    print('fetching alignments')

    algnpath = os.path.join(dtpath, algnfile)
    alignments = json.load(open(algnpath))

    if process:
        mdfile = os.path.join(dtpath, mdprefix + mode + '.pkl')
        alignments = process_alignments(alignments, mdfile, ids_to_drop)

    print('making dataset')

    data = make_XY.make_xy2(alignments)

    deepL = train_deepARG.main(data)

    # this is used to store the parameters of the neural network.
    cPickle.dump(deepL['parameters'], open(
        dtpath+"/model/"+version+"/metadata_"+mode+".pkl", "w"))
    deepL['clf'].save_params_to(dtpath+"/model/"+version+"/model_"+mode+".pkl")

if __name__ == '__main__':
    train()
