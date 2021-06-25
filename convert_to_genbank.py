import argparse
import logging
import sys

import rdflib

import sbol2
import sbol3
import tyto


FILENAME = 'iGEM_round_1_order_constructs_full_sequences.json'
OUTPUT_SBOL = 'iGEM_round_1_order_constructs_full_sequences_sbol2.xml'
OUTPUT_GENBANK = 'iGEM_round_1_order_constructs_full_sequences.gb'
OUTPUT_FASTA = 'iGEM_round_1_order_constructs_full_sequences.fasta'

###############################################################
# Converter routines

def convert_toplevel_and_dependencies(target, t):
    # If already converted, return
    converted = target.find(t.identity+'/1')
    if converted:
        return converted
    # Otherwise, convert the object
    print("Converting " + t.identity)
    if isinstance(t, sbol3.Collection):
        tnew = convert_collection_and_dependencies(target, t)
    elif isinstance(t, sbol3.Component):
        tnew = convert_component_and_dependencies(target, t)
    elif isinstance(t, sbol3.Sequence):
        tnew = convert_sequence(target, t)
    else:
        raise ValueError("Not set up to convert "+str(t))
    # Copy the shared Identified fields
    tnew.displayId = t.display_id
    tnew.name = t.name
    tnew.description = t.description
    return tnew

def convert_sequence(target, s):
    snew = sbol2.Sequence(s.identity)
    #target.add(snew) # pysbol2 bug - gets added when the sequence is set
    snew.encoding = convert_type(s.encoding)
    snew.elements = s.elements
    return snew

def convert_collection_and_dependencies(target, c):
    cnew = sbol2.Collection(c.identity)
    target.add(cnew)
    for m in c.members:
        cnew.members.append(convert_toplevel_and_dependencies(target, m.lookup()))
    return cnew

remapped_types = {
    sbol3.SBO_DNA:sbol2.BIOPAX_DNA,
    sbol3.SBO_RNA:sbol2.BIOPAX_RNA,
    sbol3.SBO_PROTEIN:sbol2.BIOPAX_PROTEIN,
    'https://identifiers.org/edam:format_1207':sbol2.SBOL_ENCODING_IUPAC #   ### BUG: pySBOL #185
}
def convert_type(type):
    return (remapped_types[type] if type in remapped_types.keys() else type)

def convert_location(location):
    if isinstance(location, sbol3.Range):
        return sbol2.Range(location.display_id, start=location.start, end=location.end)
    elif isinstance(location, sbol3.Cut):
        return sbol2.Cut(location.display_id, at=location.at)
    else:
        raise ValueError("Not set up to convert " + str(t))


def convert_component_and_dependencies(target, c):
    cnew = sbol2.ComponentDefinition(c.identity)
    target.add(cnew)
    cnew.types = [convert_type(t) for t in c.types]
    cnew.roles = c.roles
    for f in c.features:
        if isinstance(f,sbol3.SubComponent):
            subc = convert_toplevel_and_dependencies(target, f.instance_of.lookup())
            subcnew = sbol2.Component(f.display_id)
            subcnew.definition = subc.identity
            cnew.components.add(subcnew)
            if len(f.locations)>0:
                sa = sbol2.SequenceAnnotation(subcnew.displayId+"_SA")
                cnew.sequenceAnnotations.add(sa)
                sa.component = subcnew
                for l in f.locations:
                    sa.locations.add(convert_location(l))
        else:
            raise ValueError("Not set up to convert " + str(f))
    for s in c.sequences:
        snew = convert_toplevel_and_dependencies(target, s.lookup())
        cnew.sequence = snew
    # Just going to ignore the rest of the material for now
    return cnew

# kludge: it's a plasmid if either it or its direct subcomponent is
def component_is_circular_plasmid_with_inserts(component):
    plasmid = tyto.SO.get_uri_by_term('plasmid')
    return isinstance(component, sbol3.Component) and \
           any(f for f in component.features if (isinstance(f,sbol3.SubComponent) and plasmid in f.instance_of.lookup().roles))


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', metavar='SBOL3_INPUT')
    parser.add_argument('-d', '--debug', action='store_true')
    args = parser.parse_args(args)
    return args


def init_logging(debug=False):
    msg_format = '%(asctime)s %(levelname)s %(message)s'
    date_format = '%m/%d/%Y %H:%M:%S'
    level = logging.INFO
    if debug:
        level = logging.DEBUG
    logging.basicConfig(format=msg_format, datefmt=date_format, level=level)


def main(argv=None):
    args = parse_args(argv)
    init_logging(args.debug)


    ###############################################################
    # Load the SBOL3 document

    print('Reading SBOL3 document')
    input_doc = sbol3.Document()
    input_doc.read(args.infile, rdflib.util.guess_format(args.infile))

    ###############################################################
    # Add plasmid information to the description of all inserts
    # add plasmid information for all of the non-serializable

    print("Adding backbone information to insert constructs")

    plasmid = tyto.SO.get_uri_by_term('plasmid')

    plasmids = {c for c in input_doc.objects if component_is_circular_plasmid_with_inserts(c)}
    for p in plasmids:
        backbone = {f for f in p.features if (isinstance(f,sbol3.SubComponent) and plasmid in f.instance_of.lookup().roles)}
        assert len(backbone)==1, 'Only know how to convert plasmid constructs with precisely one plasmid feature: '+p.identity
        backbone = backbone.pop() # get the backbone from the set
        for c in (f for f in p.features if isinstance(f,sbol3.SubComponent)):
            #insert_description = "; Vector: "+backbone.instance_of.lookup().display_id
            #c.instance_of.lookup().description = c.instance_of.lookup().display_id + insert_description + ("; "+c.instance_of.lookup().description if c.instance_of.lookup().description else '')
            # Looks like for setting up Twist orders we'll need to use just the UID, since names have to be <32 characters on Twist; kludging for now
            c.instance_of.lookup().description = c.instance_of.lookup().display_id



    ###############################################################
    # Run the converter on all top-levels in the SBOL2 document

    print('Converting to SBOL 2')
    output_doc = sbol2.Document()
    sbol2.setHomespace('')
    sbol2.Config.setOption('sbol_typed_uris', False)
    serializable = (c for c in input_doc.objects if isinstance(c,sbol3.Component) and c.sequences)
    for toplevel in serializable:
        convert_toplevel_and_dependencies(output_doc, toplevel)

    report = output_doc.validate()
    print(report)

    output_doc.write(OUTPUT_SBOL)
    print('SBOL2 file written')

    output_doc.exportToFormat('GenBank',OUTPUT_GENBANK)
    print('GenBank file written')


if __name__ == '__main__':
    main()
