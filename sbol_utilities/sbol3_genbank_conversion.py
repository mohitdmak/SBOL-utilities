import os
import csv
import math
import sbol3
import logging
from collections import OrderedDict
from typing import Dict, List, Sequence, Union, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, Reference, \
    CompoundLocation, BeforePosition, ExactPosition, AfterPosition



class GenBank_SBOL3_Converter:
    """Main Converter class handling offline, direction conversion of files SBOL3 files to and from GenBank files"""
    # dictionaries to store feature lookups for terms in GenBank and SO ontologies
    gb2so_map = {}
    so2gb_map = {}
    ## Conversion Constants :
    # TODO: Temporarily assuming only dna components to be dealt with in genbank files
    COMP_TYPES = [sbol3.SBO_DNA]
    # TODO: Temporarily assuming components to only have the engineered_region role
    COMP_ROLES = [sbol3.SO_ENGINEERED_REGION]
    # TODO: Temporarily encoding sequnce objects in IUPUC mode only
    SEQUENCE_ENCODING = sbol3.IUPAC_DNA_ENCODING
    # BIO_STRAND constants, which server as the GenBank counterparts to SBOL3's inline and reverse orientations
    BIO_STRAND_FORWARD = 1
    BIO_STRAND_REVERSE = -1
    # Mapping int types to types of locationPositions in GenBank (Before/After/Exact)
    SBOL_LOCATION_POSITION = {BeforePosition: 0, ExactPosition: 1, AfterPosition: 2}
    GENBANK_LOCATION_POSITION = {0: BeforePosition, 1: ExactPosition, 2: AfterPosition}
    # Default value for the "sequence_version" annotation in GenBank files
    DEFAULT_GB_SEQ_VERSION = 1
    # Default terms for SBOL3 and GenBank in case the feature lookup from respective dictionaries does not yield any ontology term
    DEFAULT_SO_TERM = "SO:0000110"
    DEFAULT_GB_TERM = "misc_feature"
    # Namespace to be used be default if not provided, and also for all unit tests related to this converter
    TEST_NAMESPACE = "https://test.sbol3.genbank/"
    # File locations for required CSV data files which store the ontology term translations between GenBank and SO ontologies
    GB2SO_MAPPINGS_CSV = os.path.join(os.path.dirname(os.path.realpath(__file__)), "gb2so.csv")
    SO2GB_MAPPINGS_CSV = os.path.join(os.path.dirname(os.path.realpath(__file__)), "so2gb.csv")


    def __init__(self) -> None:
        """While instantiating an instance of the converter, required builders
        must be registered in order to accurately parse modified or new SBOL3 class objects
        """
        def build_component_genbank_extension(*, identity, type_uri) -> GenBank_SBOL3_Converter.Component_GenBank_Extension:
            """A builder function to be called by the SBOL3 parser
            when it encounters a Component in an SBOL file.
            :param identity: identity for new component class instance to have
            :param type_uri: type_uri for new component class instance to have
            """
            # `types` is required and not known at build time.
            # Supply a missing value to the constructor, then clear
            # the missing value before returning the built object.
            obj = self.Component_GenBank_Extension(identity=identity, types=[sbol3.PYSBOL3_MISSING], type_uri=type_uri)
            # Remove the placeholder value
            obj.clear_property(sbol3.SBOL_TYPE)
            return obj

        def build_feature_qualifiers_extension(*, identity, type_uri) -> GenBank_SBOL3_Converter.Feature_GenBank_Extension:
            """A builder function to be called by the SBOL3 parser
            when it encounters a SequenceFeature in an SBOL file.
            :param identity: identity for new sequence feature class instance to have
            :param type_uri: type_uri for new sequence feature class instance to have
            """
            # `types` is required and not known at build time.
            # Supply a missing value to the constructor, then clear
            # the missing value before returning the built object.
            obj = self.Feature_GenBank_Extension(identity=identity, type_uri=type_uri)
            # Remove the placeholder value
            obj.clear_property(sbol3.SBOL_TYPE)
            return obj

        def build_range_extension(*, identity, type_uri) -> GenBank_SBOL3_Converter.Range_GenBank_Extension:
            """A builder function to be called by the SBOL3 parser
            when it encounters a Range location in an SBOL file.
            :param identity: identity for new sequence feature class instance to have
            :param type_uri: type_uri for new sequence feature class instance to have
            """
            # `types` is required and not known at build time.
            # Supply a missing value to the constructor, then clear
            # the missing value before returning the built object.
            obj = self.Range_GenBank_Extension(identity=identity, type_uri=type_uri)
            # Remove the placeholder value
            obj.clear_property(sbol3.SBOL_TYPE)
            return obj

        def build_custom_reference_property(*, identity, type_uri) -> GenBank_SBOL3_Converter.CustomReferenceProperty:
            """A builder function to be called by the SBOL3 parser
            when it encounters a CustomReferenceProperty Toplevel object in an SBOL file.
            :param identity: identity for custom reference property instance to have
            :param type_uri: type_uri for custom reference property instance to have
            """
            obj = self.CustomReferenceProperty(identity=identity, type_uri=type_uri)
            return obj

        def build_custom_structured_comment_property(*, identity, type_uri) -> GenBank_SBOL3_Converter.CustomStructuredCommentProperty:
            """A builder function to be called by the SBOL3 parser
            when it encounters a CustomStructuredCommentProperty Toplevel object in an SBOL file.
            :param identity: identity for custom comment property instance to have
            :param type_uri: type_uri for custom comment property instance to have
            """
            obj = self.CustomStructuredCommentProperty(identity=identity, type_uri=type_uri)
            return obj

        # Register the builder function so it can be invoked by
        # the SBOL3 parser to build objects with a Component type URI
        sbol3.Document.register_builder(sbol3.SBOL_COMPONENT, build_component_genbank_extension)
        # Register the buildre function for custom reference properties
        sbol3.Document.register_builder(self.CustomReferenceProperty.CUSTOM_REFERENCE_NS, build_custom_reference_property)
        # Register the buildre function for custom structured comment properties
        sbol3.Document.register_builder(self.CustomStructuredCommentProperty.CUSTOM_STRUCTURED_COMMENT_NS, build_custom_structured_comment_property)
        # # Register the builder function so it can be invoked by
        # # the SBOL3 parser to build objects with a SequenceFeature type URI
        sbol3.Document.register_builder(sbol3.SBOL_SEQUENCE_FEATURE, build_feature_qualifiers_extension)
        # # Register the builder function so it can be invoked by
        # # the SBOL3 parser to build objects with a Range type URI
        sbol3.Document.register_builder(sbol3.SBOL_RANGE, build_range_extension)


    class CustomReferenceProperty(sbol3.CustomTopLevel):
        """Serves to store information and annotations for 'Reference' objects in 
        GenBank file to SBOL3 while parsing so that it may be retrieved back in a round trip
        :extends: sbol3.CustomTopLevel class
        """
        CUSTOM_REFERENCE_NS = "http://www.ncbi.nlm.nih.gov/genbank#reference"
        def __init__(self, type_uri=CUSTOM_REFERENCE_NS, identity="customReferenceProperty"):
            super().__init__(identity, type_uri)
            self.authors    = sbol3.TextProperty(self, f"{self.CUSTOM_REFERENCE_NS}#authors"   , 0, 1)
            self.comment    = sbol3.TextProperty(self, f"{self.CUSTOM_REFERENCE_NS}#comment"   , 0, 1)
            self.journal    = sbol3.TextProperty(self, f"{self.CUSTOM_REFERENCE_NS}#journal"   , 0, 1)
            self.consrtm    = sbol3.TextProperty(self, f"{self.CUSTOM_REFERENCE_NS}#consrtm"   , 0, 1)
            self.title      = sbol3.TextProperty(self, f"{self.CUSTOM_REFERENCE_NS}#title"     , 0, 1)
            self.medline_id = sbol3.TextProperty(self, f"{self.CUSTOM_REFERENCE_NS}#medline_id", 0, 1)
            self.pubmed_id  = sbol3.TextProperty(self, f"{self.CUSTOM_REFERENCE_NS}#pubmed_id" , 0, 1)
            # stores the display id of parent component for a particular CustomReferenceProperty object
            self.component  = sbol3.TextProperty(self, f"{self.CUSTOM_REFERENCE_NS}#component" , 0, 1)
            # TODO: support cut locations?
            # there can be multiple locations described for a reference, thus upper bound needs to be > 1 in order to use ListProperty
            self.location = sbol3.OwnedObject(self, f"{self.CUSTOM_REFERENCE_NS}#location", 0, math.inf, type_constraint=sbol3.Range)


    class CustomStructuredCommentProperty(sbol3.CustomTopLevel):
        """Serves to store information and annotations for 'Structured_Comment' objects in 
        GenBank file to SBOL3 while parsing so that it may be retrieved back in a round trip
        Complete reference available at: https://www.ncbi.nlm.nih.gov/genbank/structuredcomment/
        :extends: sbol3.CustomTopLevel class
        """
        CUSTOM_STRUCTURED_COMMENT_NS = "http://www.ncbi.nlm.nih.gov/genbank#structured_comment"
        def __init__(self, type_uri=CUSTOM_STRUCTURED_COMMENT_NS, identity="customStructuredCommentProperty"):
            super().__init__(identity, type_uri)
            self.heading    = sbol3.TextProperty(self, f"{self.CUSTOM_STRUCTURED_COMMENT_NS}#heading"   , 0, 1)
            # stores the display id of parent component for a particular CustomReferenceProperty object
            self.component  = sbol3.TextProperty(self, f"{self.CUSTOM_STRUCTURED_COMMENT_NS}#component" , 0, 1)
            # there can be multiple key/values described for a structured_comment, thus upper bound needs to be > 1 in order to use ListProperty
            self.structured_keys   = sbol3.TextProperty(self, f"{self.CUSTOM_STRUCTURED_COMMENT_NS}#structuredKeys", 0, math.inf)
            self.structured_values = sbol3.TextProperty(self, f"{self.CUSTOM_STRUCTURED_COMMENT_NS}#structuredValues", 0, math.inf)


    class Feature_GenBank_Extension(sbol3.SequenceFeature):
        """Overrides the sbol3 SequenceFeature class to include fields to directly read and write 
        qualifiers of GenBank features not storeable in any SBOL3 datafield.
        :extends: sbol3.SequenceFeature class
        """
        GENBANK_FEATURE_QUALIFIER_NS = "http://www.ncbi.nlm.nih.gov/genbank#featureQualifier"
        def __init__(self, locations: List[sbol3.Location] = [], **kwargs) -> None:
            # instantiating sbol3 SequenceFeature object
            super().__init__(locations=locations, **kwargs)
            # Setting properties for GenBank's qualifiers not settable in any SBOL3 field.
            self.qualifier_key      = sbol3.TextProperty(self, f"{self.GENBANK_FEATURE_QUALIFIER_NS}#key"  , 0, math.inf)
            self.qualifier_value    = sbol3.TextProperty(self, f"{self.GENBANK_FEATURE_QUALIFIER_NS}#value", 0, math.inf)


    class Range_GenBank_Extension(sbol3.Range):
        """Overrides the sbol3 Range class to include fields to store the  
        start and end position types (AfterPostion / BeforePosition / ExactPosition).
        :extends: sbol3.Range class
        """
        GENBANK_RANGE_NS = "http://www.ncbi.nlm.nih.gov/genbank#locationPosition"
        def __init__(self, sequence: sbol3.Sequence = sbol3.Sequence("autoCreatedSequence"), start: int = 0, end: int = 0, **kwargs) -> None:
            # instantiating sbol3 SequenceFeature object
            super().__init__(sequence = sequence, start = start, end = end, **kwargs)
            # Setting properties for GenBank's location position not settable in any SBOL3 field.
            self.start_position = sbol3.IntProperty(self, f"{self.GENBANK_RANGE_NS}#start", 0, 1)
            self.end_position   = sbol3.IntProperty(self, f"{self.GENBANK_RANGE_NS}#end"  , 0, 1)


    class Component_GenBank_Extension(sbol3.Component):
        """Overrides the sbol3 Component class to include fields to directly read and write 
        extraneous properties of GenBank not storeable in any SBOL3 datafield.
        :extends: sbol3.Component class
        """
        GENBANK_EXTRA_PROPERTY_NS = "http://www.ncbi.nlm.nih.gov/genbank"
        def __init__(self, identity: str, types: Optional[Union[str, Sequence[str]]], **kwargs) -> None:
            # instantiating sbol3 component object
            super().__init__(identity=identity, types=types, **kwargs)
            # Setting properties for GenBank's extraneous properties not settable in any SBOL3 field.
            self.genbank_seq_version   = sbol3.IntProperty(self,  f"{self.GENBANK_EXTRA_PROPERTY_NS}#seq_version", 0, 1)
            self.genbank_date          = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#date"       , 0, 1)
            self.genbank_division      = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#division"   , 0, 1)
            self.genbank_locus         = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#locus"      , 0, 1)
            self.genbank_molecule_type = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#molecule"   , 0, 1)
            self.genbank_organism      = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#organism"   , 0, 1)
            self.genbank_source        = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#source"     , 0, 1)
            self.genbank_topology      = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#topology"   , 0, 1)
            self.genbank_gi            = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#gi"         , 0, 1)
            self.genbank_comment       = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#comment"    , 0, 1)
            self.genbank_dblink        = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#dbxrefs"    , 0, 1)
            self.genbank_record_id     = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#id"         , 0, 1)
            # TODO : add note linking issue here
            self.genbank_taxonomy      = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#taxonomy"   , 0, 1)
            self.genbank_keywords      = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#keywords"   , 0, 1)
            # there can be multiple accessions, thus upper bound needs to be > 1 in order to use TextListProperty
            self.genbank_accessions    = sbol3.TextProperty(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#accessions", 0, math.inf)
            self.fuzzy_features        = sbol3.OwnedObject(self, f"{self.GENBANK_EXTRA_PROPERTY_NS}#fuzzyFeatures", 0, math.inf, type_constraint=sbol3.SequenceFeature)


    def create_GB_SO_role_mappings(self, gb2so_csv: str = GB2SO_MAPPINGS_CSV, so2gb_csv: str = SO2GB_MAPPINGS_CSV,
                                   convert_gb2so: bool = True, convert_so2gb: bool = True) -> int:
        """Reads 2 CSV Files containing mappings for converting between GenBank and SO ontologies roles
        :param gb2so_csv: path to read genbank to so conversion csv file
        :param so2gb_csv: path to read so to genbank conversion csv file
        :param convert_gb2so: bool stating whether to read csv for gb2so mappings
        :param convert_so2gb: bool stating whether to read csv for so2gb mappings
        :return: int 1 / 0 denoting the status of whether the mappings were created and stored in dictionaries
        """
        if convert_gb2so:
            logging.info(f"Parsing {gb2so_csv} for GenBank to SO ontology mappings.")
            try:
                with open(gb2so_csv, mode="r") as csv_file:
                    csv_reader = csv.DictReader(csv_file)
                    for row in csv_reader:
                        self.gb2so_map[row["GenBank_Ontology"]] = row["SO_Ontology"]
            except FileNotFoundError:
                logging.error(f"No GenBank to SO Ontology Mapping CSV File Exists!")
                return 0
        if convert_so2gb:
            logging.info(f"Parsing {so2gb_csv} for SO to GenBank ontology mappings.")
            try:
                with open(so2gb_csv, mode="r") as csv_file:
                    csv_reader = csv.DictReader(csv_file)
                    for row in csv_reader:
                        self.so2gb_map[row["SO_Ontology"]] = row["GenBank_Ontology"]
            except FileNotFoundError:
                logging.error(f"No SO to Genbank Ontology Mapping CSV File Exists!")
                return 0
        return 1


    def convert_genbank_to_sbol3(self, gb_file: str, sbol3_file: str = "sbol3.nt", namespace: str = TEST_NAMESPACE,
                                 write: bool = False) -> sbol3.Document:
        """Convert a GenBank document on disk into an SBOL3 document
        The GenBank document is parsed using BioPython, and corresponding objects of SBOL3 document are created

        :param gb_file: path to read GenBank file from
        :param sbol3_file: path to write SBOL3 file to, if write set to true
        :param namespace: URIs of Components will be set to {namespace}/{genbank_id},
                          defaults to "https://test.sbol3.genbank/"
        :param write: writes the generated sbol3 document in SORTED_NTRIPLES
                          format to provided sbol3_file path
        :return: SBOL3 document containing converted materials
        """
        # create sbol3 document, and record parser handler for gb file
        sbol3.set_namespace(namespace)
        doc = sbol3.Document()
        # create updated py dict to store mappings between gb and so ontologies
        logging.info(
            "Creating GenBank and SO ontologies mappings for sequence feature roles"
        )
        map_created = self.create_GB_SO_role_mappings(
            gb2so_csv=self.GB2SO_MAPPINGS_CSV, convert_so2gb=False
        )
        if not map_created:
            # TODO: Need better SBOL3-GenBank specific error classes in future
            raise ValueError(
                "Required CSV data files are not present in your package.\n    Please reinstall the sbol_utilities package.\n \
                Stopping current conversion process.\n    Reverting to legacy converter if new Conversion process is not forced."
            )
        # access records by parsing gb file using SeqIO class
        logging.info(
            f"Parsing Genbank records using SeqIO class.\n    Using GenBank file {gb_file}"
        )
        for record in list(SeqIO.parse(gb_file, "genbank").records):
            # TODO: Currently we assume only linear or circular topology is possible
            logging.info(f"Parsing record - `{record.id}` in genbank file.")
            topology = "linear"
            if "topology" in record.annotations:
                topology = record.annotations["topology"]
            # sometimes topology is specified in the 'data_file_division' field
            elif record.annotations['data_file_division'] in ['circular', 'linear']:
                topology = record.annotations['data_file_division']
            if topology == "linear":
                extra_comp_types = [sbol3.SO_LINEAR]
            else:
                extra_comp_types = [sbol3.SO_CIRCULAR]
            # creating component extended Component class to include GenBank extraneous properties
            comp = self.Component_GenBank_Extension(
                identity=record.name,
                types=self.COMP_TYPES + extra_comp_types,
                roles=self.COMP_ROLES,
                description=record.description,
            )
            doc.add(comp)

            # TODO: Currently we use a fixed method of encoding (IUPAC)
            seq = sbol3.Sequence(
                identity=record.name + "_sequence",
                elements=str(record.seq.lower()),
                encoding=self.SEQUENCE_ENCODING,
            )
            doc.add(seq)
            comp.sequences = [seq]

            # Setting properties for GenBank's extraneous properties not settable in any SBOL3 field.
            self._store_extra_properties_in_sbol3(doc, comp, seq, record)

            # create all sequence features, and tag all encountered feature qualifiers via extended Feature_GenBank_Extension class
            self._handle_features_gb_to_sbol(record, comp, seq)

        if write:
            logging.info(
                f"Writing created sbol3 document to disk in sorted ntriples format.\n    With path {sbol3_file}"
            )
            doc.write(fpath=sbol3_file, file_format=sbol3.SORTED_NTRIPLES)
        return doc


    def convert_sbol3_to_genbank(self, sbol3_file: str, doc: sbol3.Document = None, gb_file: str = "genbank.out",
                                 # write: bool = False) -> List[SeqRecord]:
                                 write: bool = False) -> Dict:
        """Convert a SBOL3 document on disk into a GenBank document
        The GenBank document is made using an array of SeqRecords using BioPython, by parsing SBOL3 objects

        :param sbol3_file: path to read SBOL3 file from
        :param gb_file: path to write GenBank file to, if write set to true
        :param write: writes the generated genbank document to provided path
        :return: Array of SeqRecord objects which comprise the generated GenBank document
        """
        if not doc:
            doc = sbol3.Document()
            doc.read(sbol3_file)
        seq_records = []
        logging.info(
            "Creating GenBank and SO ontologies mappings for sequence feature roles"
        )
        # create logs dict to be returned as conversion status of the SBOL3 file provided
        logs: Dict[sbol3.TopLevel, bool] = {}
        # create dict to link component with their respective Reference property objects
        references: Dict[sbol3.Component, List[sbol3.CustomTopLevel]] = {}
        for obj in doc.objects:  
            if isinstance(obj, sbol3.CustomTopLevel) and obj.type_uri == self.CustomReferenceProperty.CUSTOM_REFERENCE_NS:
                component_object = doc.find(str(obj.component))
                if component_object and isinstance(component_object, sbol3.Component):
                    references[component_object] = [obj] if component_object not in references else references[component_object] + [obj]
                # TODO: Raise error here
                # else:
        # create dict to link component with their respective structured comment objects
        structured_comments: Dict[sbol3.Component, List[sbol3.CustomTopLevel]] = {}
        for obj in doc.objects:  
            if isinstance(obj, sbol3.CustomTopLevel) and obj.type_uri == self.CustomStructuredCommentProperty.CUSTOM_STRUCTURED_COMMENT_NS:
                component_object = doc.find(str(obj.component))
                if component_object and isinstance(component_object, sbol3.Component):
                    structured_comments[component_object] = [obj] if component_object not in structured_comments else structured_comments[component_object] + [obj]
                # TODO: Raise error here
                # else:

        # create updated py dict to store mappings between gb and so ontologies
        map_created = self.create_GB_SO_role_mappings(
            so2gb_csv=self.SO2GB_MAPPINGS_CSV, convert_gb2so=False
        )
        if not map_created:
            # TODO: Need better SBOL3-GenBank specific error classes in future
            raise ValueError(
                f"Required CSV data files are not present in your package.\n    Please reinstall the sbol_utilities package.\n \
                Stopping current conversion process.\n    Reverting to legacy converter if new Conversion process is not forced."
            )
        # consider sbol3 objects which are components
        logging.info(f"Parsing SBOL3 Document components using SBOL3 Document: \n{doc}")
        for obj in doc.objects:
            if isinstance(obj, sbol3.TopLevel):
                # create a key for the top level object if it is not already parsed
                if obj not in logs:
                    logs[obj] = False
            if isinstance(obj, sbol3.Component):
                logging.info(f"Parsing component - `{obj.display_id}` in sbol3 document.")
                # NOTE: A single component/record cannot have multiple sequences
                seq = None # If no sequence is found for a component
                if obj.sequences and len(obj.sequences) == 1:
                    if doc.find(obj.sequences[0]):
                        obj_seq = doc.find(obj.sequences[0])
                        seq = Seq(obj_seq.elements.upper())
                        # mark the status of this top level sequence object as parsed and converted
                        if isinstance(obj_seq, sbol3.TopLevel): 
                            logs[obj_seq] = True
                elif len(obj.sequences) > 1:
                    raise ValueError(f"Component `{obj.display_id}` of given SBOL3 document has more than 1 sequnces \n \
                    (`{len(obj.sequences)}`). This is invalid; a component may only have 1 or 0 sequences.")
                seq_rec = SeqRecord(
                    seq=seq,
                    id=obj.display_id,
                    description=obj.description,
                    name=obj.display_id,
                )
                # Resetting extraneous genbank properties from extended component-genbank class
                self._reset_extra_properties_in_genbank(obj, seq_rec, references, structured_comments)

                seq_rec_features = []
                # recreate all sequence features, and tag all encountered feature qualifiers via extended Feature_GenBank_Extension class
                self._handle_features_sbol_to_gb(seq_rec, obj, seq_rec_features)

                # mark the top level component object as parsed and converter
                logs[obj] = True
                seq_records.append(seq_rec)
        # writing generated genbank document to disk at path provided
        if write:
            logging.info(
                f"Writing created genbank file to disk.\n    With path {gb_file}"
            )
            SeqIO.write(seq_records, gb_file, "genbank")
        return {"status": logs, "seqrecords": seq_records}


    def _store_extra_properties_in_sbol3(self, doc: sbol3.Document, comp: Component_GenBank_Extension, 
                                       seq: sbol3.Sequence, record: SeqRecord) -> None:
        """Helper function for setting properties for GenBank's extraneous properties not directly settable in any SBOL3 field,
        by using a modified, extended SBOL3 Component class, and a new CustomReferenceProperty TopLevel class.
        :param doc: The sbol3 document to store the contents in
        :param comp: Instance of the extended SBOL3 Component class (Component_GenBank_Extension)
        :param seq: The Sequence used in the GenBank record corresponding to sbol3 comp
        :param record: GenBank SeqRecord instance for the record which contains extra properties
        """
        comp.genbank_record_id = record.id
        # set dblinks from the dbxrefs property of biopython
        if record.dbxrefs:
            # dbxrefs are parsed in a list by biopython from `record.dbxrefs`; we are storing them as a flat string
            # so as to maintain order. Thus, we are creating a custom delimetere of `::`, by which we shall seperate
            # individual dbxrefs in the string and later split them to a list while resetting them in genbank
            comp.genbank_dblink = "::".join(record.dbxrefs)
        for annotation in record.annotations:
            # Sending out warnings for genbank info not storeable in sbol3
            logging.warning(
                f"Extraneous information not directly storeable in SBOL3 - {annotation}: {record.annotations[annotation]}"
            )
            # 1. GenBank Record Date
            if annotation == 'date':
                comp.genbank_date = record.annotations['date']
            # 2. GenBank Record Division
            elif annotation == 'data_file_division':
                # FIX for iGEM files not having data file division but topology stored in its key
                if record.annotations['data_file_division'] in ['circular', 'linear']:
                    comp.genbank_topology = record.annotations['data_file_division']
                else: comp.genbank_division = record.annotations['data_file_division']
            # 3. GenBank Record Keywords
            elif annotation == 'keywords':
                comp.genbank_keywords = ",".join(record.annotations['keywords'])
            # 4. GenBank Record Molecule Type
            elif annotation == 'molecule_type':
                comp.genbank_molecule_type = record.annotations['molecule_type']
            # 5. GenBank Record Organism
            elif annotation == 'organism':
                comp.genbank_organism = record.annotations['organism']
            # 6. GenBank Record Source
            elif annotation == 'source':
                comp.genbank_source = record.annotations['source']
            # 7. GenBank Record Taxonomy
            elif annotation == 'taxonomy':
                comp.genbank_taxonomy = ",".join(record.annotations['taxonomy'])
            # 8. GenBank Record Topology
            elif annotation == 'topology':
                comp.genbank_topology = record.annotations['topology']
            # 9. GenBank Record GI Property
            elif annotation == 'gi':
                comp.genbank_gi = record.annotations['gi']
            # 10. GenBank Record Accessions
            elif annotation == 'accessions':
                comp.genbank_accessions = sorted(record.annotations['accessions'])
            # 11. GenBank Sequence Version
            elif annotation == 'sequence_version':
                comp.genbank_seq_version = record.annotations['sequence_version']
            # 12. GenBank Record References
            elif annotation == 'references':
                for ind, reference in enumerate(record.annotations['references']):
                    # create a custom reference property instance for each reference
                    custom_reference = self.CustomReferenceProperty(identity = comp.identity + f"/Reference_{ind}")
                    custom_reference.authors = reference.authors
                    custom_reference.comment = reference.comment
                    custom_reference.journal = reference.journal
                    custom_reference.title = reference.title
                    custom_reference.consrtm = reference.consrtm
                    custom_reference.medline_id = reference.medline_id
                    custom_reference.pubmed_id = reference.pubmed_id
                    for gb_loc in reference.location:
                        feat_loc_orientation = sbol3.SO_FORWARD
                        if gb_loc.strand == -1:
                            feat_loc_orientation = sbol3.SO_REVERSE
                        if gb_loc.start == gb_loc.end:
                            locs = sbol3.Cut(sequence=seq, at=int(gb_loc.start), orientation=feat_loc_orientation)
                        else:
                            locs = sbol3.Range(sequence=seq, start=int(gb_loc.start), end=int(gb_loc.end), orientation=feat_loc_orientation)
                        custom_reference.location.append(locs)
                    # link the parent component for each custom reference property objects
                    if comp.display_id:
                        custom_reference.component = comp.display_id
                    # TODO: Raise error, no name for component
                    # else:
                    doc.add(custom_reference)
            # 14. GenBank Record Comment
            elif annotation == 'comment':
                comp.genbank_comment = record.annotations['comment']
            # 15. GenBank Record Structured comments
            elif annotation == 'structured_comment':
                identity_ind = 1
                for heading in record.annotations['structured_comment']:
                    structured_comment_object = self.CustomStructuredCommentProperty(identity = comp.identity + f"/StructuredComment_{identity_ind}") 
                    identity_ind += 1
                    if comp.display_id:
                        structured_comment_object.component = comp.display_id
                    structured_comment_object.heading = heading
                    structured_dict = record.annotations['structured_comment'][heading]
                    key_value_ind = 1
                    for key in structured_dict:
                        # NOTE: if storing list as string for keys and values both, have a check of them having same length user uses our delimeter while writing
                        structured_comment_object.structured_keys.append(f"{key_value_ind}::{key}")
                        structured_comment_object.structured_values.append(f"{key_value_ind}::{structured_dict[key]}")
                        key_value_ind += 1
                    doc.add(structured_comment_object)
            else:
                raise ValueError(f"The annotation `{annotation}` in the GenBank record `{record.id}`\n \
                                    is not recognized as a standard annotation.")
        # TODO: BioPython's parsing doesn't explicitly place a "locus" datafield?
        # 13. GenBank Record Locus
        comp.genbank_locus = record.name


    def _reset_extra_properties_in_genbank(self, obj: sbol3.Component, seq_rec: SeqRecord, references: Dict[sbol3.Component, List[sbol3.CustomTopLevel]], structured_comments: Dict[sbol3.Component, List[sbol3.CustomTopLevel]]) -> None:
        """Helper function for resetting properties for GenBank's extraneous properties from SBOL3 Document's properties,
        by using a modified, extended SBOL3 Component class, and a new CustomReferenceProperty TopLevel class.
        :param obj: SBOL3 component, extra properties would be stored within it if its an instance of the extended SBOL3 Component class
        :param seq_rec: GenBank SeqRecord instance for the record in which to reset extra properties
        :param references: Dictionary linking SBOL3 components to their respective list of CustomReferenceProperty objects
        """
        if isinstance(obj, self.Component_GenBank_Extension):
            seq_rec.id = obj.genbank_record_id
            # set db links using dbxrefs property of biopython
            if obj.genbank_dblink:
                # NOTE: see comment on `_store_extra_properties_in_sbol3`'s dbxrefs section, wherein we describe how '::' 
                # is used as a delimeter to store the dbxrefs list as a string to maintain order. Here, we split the string 
                # by the same delimeter to restore the list in resetting GenBank properties.
                seq_rec.dbxrefs = str(obj.genbank_dblink).split("::")
            # 1. GenBank Record Date
            seq_rec.annotations['date'] = obj.genbank_date
            # 2. GenBank Record Division
            seq_rec.annotations['data_file_division'] = obj.genbank_division
            # 3. GenBank Record Keywords
            # seq_rec.annotations['keywords'] = sorted(list(obj.genbank_keywords))
            if obj.genbank_keywords: seq_rec.annotations['keywords'] = str(obj.genbank_keywords).split(",")
            # 4. GenBank Record Molecule Type
            seq_rec.annotations['molecule_type'] = obj.genbank_molecule_type
            # 5. GenBank Record Organism
            seq_rec.annotations['organism'] = obj.genbank_organism
            # 6. GenBank Record Source
            # FIXME: Apparently, if a default source was used during in the GenBank file
            #        during conversion of GenBank -> SBOL, component.genbank_source is "", 
            #        and while plugging it back in during conversion of SBOL -> GenBank, it
            #        simply prints "", whereas the default "." should have been printed
            if obj.genbank_source != "": seq_rec.annotations['source'] = obj.genbank_source
            # 7. GenBank Record taxonomy
            # TODO : link gh issue for note below
            # FIXME: Even though component.genbank_taxonomy is stored in sorted order, it 
            #        becomes unsorted while retrieving from the sbol file
            if obj.genbank_taxonomy: seq_rec.annotations['taxonomy'] = str(obj.genbank_taxonomy).split(",")
            # 8. GenBank Record Topology
            seq_rec.annotations['topology'] = obj.genbank_topology
            # 9. GenBank Record GI Property
            if obj.genbank_gi: seq_rec.annotations['gi'] = obj.genbank_gi
            # 10. GenBank Record Accessions
            seq_rec.annotations['accessions'] = sorted(list(obj.genbank_accessions))
            # 11. GenBank Sequence Version
            seq_rec.annotations['sequnce_version'] = obj.genbank_seq_version
            # 12. GenBank Record References
            if obj in references:
                # if sbol3 object has references
                record_references = []
                for reference in references[obj]:
                    reference_object = Reference()
                    reference_object.title = reference.title
                    reference_object.authors = reference.authors
                    reference_object.comment = reference.comment
                    reference_object.journal = reference.journal
                    reference_object.consrtm = reference.consrtm
                    reference_object.pubmed_id = reference.pubmed_id
                    reference_object.medline_id = reference.medline_id
                    for obj_feat_loc in reference.location:
                        feat_strand = self.BIO_STRAND_FORWARD
                        # feature strand value which denotes orientation of the location of the feature
                        # By default its 1 for SO_FORWARD orientation of sbol3 feature location, and -1 for SO_REVERSE
                        if obj_feat_loc.orientation == sbol3.SO_REVERSE:
                            feat_strand = self.BIO_STRAND_REVERSE
                        # elif obj_feat_loc.orientation != sbol3.SO_FORWARD:
                        #     raise ValueError(f"Location orientation: `{obj_feat_loc.orientation}` for feature: \n \
                        #     `{obj_feat.name}` of component: `{obj.display_id}` is not a valid orientation.\n \
                        #     Valid orientations are `{sbol3.SO_FORWARD}`, `{sbol3.SO_REVERSE}`")
                        # TODO: Raise custom converter class ERROR for `else:`
                        feat_loc_object = FeatureLocation(
                            start=obj_feat_loc.start,
                            end=obj_feat_loc.end,
                            strand=feat_strand,
                        )
                        reference_object.location.append(feat_loc_object)
                    record_references.append(reference_object)
                seq_rec.annotations['references'] = record_references
            # 14. GenBank Record Comments
            if obj.genbank_comment:
                seq_rec.annotations['comment'] = obj.genbank_comment
            # 15. GenBank Record Structured Comments
            if obj in structured_comments:
                comment_annotation = OrderedDict()
                for structured_comment in structured_comments[obj]:
                    structured_comment_object = OrderedDict()
                    total_keys = len(structured_comment.structured_keys)
                    structured_keys = sorted(list(structured_comment.structured_keys), key = lambda t: int(t.split("::", 1)[0]))
                    structured_values = sorted(list(structured_comment.structured_values), key = lambda t: int(t.split("::", 1)[0]))
                    for ind in range(total_keys):
                        key = structured_keys[ind].split("::", 1)[1]
                        value = structured_values[ind].split("::", 1)[1]
                        structured_comment_object[key] = value
                    comment_annotation[structured_comment.heading] = structured_comment_object
                seq_rec.annotations['structured_comment'] = comment_annotation
        # TODO: No explicit way to set locus via BioPython?
        # 13. GenBank Record Locus
        # TODO: temporalily hardcoding version as "1"
        seq_rec.annotations["sequence_version"] = self.DEFAULT_GB_SEQ_VERSION

    
    def _handle_features_gb_to_sbol(self, record: SeqRecord, comp: Component_GenBank_Extension, seq: sbol3.Sequence) -> None:
        """Helper function for setting sequence features and their qualifiers to SBOL,
        by using a modified, extended SBOL3 Sequence Feature class - Feature_GenBank_Extension.
        :param record: GenBank SeqRecord instance for the record which contains sequnce features
        :param comp: Instance of the SBOL3 Component 
        :param seq: The Sequence used in the GenBank record corresponding to sbol3 comp
        """
        # parse if genbank record has any features
        if not record.features:
            return
        comp.features = []
        for ind, gb_feat in enumerate(record.features):
            feat_locations = []
            fuzzy_feature = False
            feat_name = f"_converted_feature_{ind}"
            if "label" in gb_feat.qualifiers:
                feat_name = gb_feat.qualifiers["label"][0]
            logging.info(
                f"Parsing feature `{feat_name}` for record `{record.id}`"
            )
            for gb_loc in gb_feat.location.parts:
                # Default orientation is "inline" except if complement is specified via strand
                feat_loc_orientation = sbol3.SO_FORWARD
                if gb_loc.strand == -1:
                    feat_loc_orientation = sbol3.SO_REVERSE
                # create "Range/Cut" FeatureLocation by parsing genbank record location
                # Create a cut or range as featurelocation depending on whether location is specified as
                # Cut (eg: "n^n+1", parsed as [n:n] by biopython) or Range (eg: "n..m", parsed as [n:m] by biopython)
                if gb_loc.start == gb_loc.end:
                    locs = sbol3.Cut(
                        sequence=seq,
                        at=int(gb_loc.start),
                        orientation=feat_loc_orientation,
                    )
                else:
                    locs = self.Range_GenBank_Extension(
                        sequence=seq,
                        start=int(gb_loc.start),
                        end=int(gb_loc.end),
                        orientation=feat_loc_orientation,
                    )
                    # storing location types in IntProperties of SBOL3
                    locs.end_position = self.SBOL_LOCATION_POSITION[type(gb_loc.end)]
                    locs.start_position = self.SBOL_LOCATION_POSITION[type(gb_loc.start)]
                    # if any of the location endpoints of a feature (start/end) has a fuzzy end
                    # (i.e not Exact position) like BeforePosition/AfterPostion, we mark the 
                    # feature as a 'fuzzy_feature' which decides whether to store the feature or not
                    if not fuzzy_feature and locs.end_position != 1 or locs.start_position != 1:
                        fuzzy_feature = True
                feat_locations.append(locs)
            # Obtain sequence feature role from gb2so mappings
            feat_role = sbol3.SO_NS[:-3]
            if self.gb2so_map.get(gb_feat.type):
                feat_role += self.gb2so_map[gb_feat.type]
            else:
                logging.warning(f"Feature type: `{gb_feat.type}` for feature: `{gb_feat.qualifiers['label'][0]}` \n \
                of record: `{record.name}` has no corresponding ontology term for SO, using the default SO term, {self.DEFAULT_SO_TERM}")
                feat_role += self.DEFAULT_SO_TERM
            # assign feature orientation based on the strand value in genbank feature
            feat_orientation = sbol3.SO_FORWARD
            if gb_feat.strand == -1:
                feat_orientation = sbol3.SO_REVERSE
            feat = self.Feature_GenBank_Extension(
                locations=feat_locations,
                roles=[feat_role],
                # name=gb_feat.qualifiers["label"][0],
                name=feat_name,
                orientation=feat_orientation
            )
            # store qualifiers key value pairs
            for ind, qualifier in enumerate(gb_feat.qualifiers):
                feat.qualifier_key.append(f"{ind}:" + qualifier)
                feat.qualifier_value.append(f"{ind}:" + gb_feat.qualifiers[qualifier][0])
            # if feature has any fuzzy location, since SBOL does not support storing such location endpoints,
            # instead of presenting incomplete/incorrect information to users, we would store the feature assume
            # a property of the Extended GenBank Component class, instead as a feature of the component.
            # See: issue -> 
            # Once the above issue gets addressed, we can remove the 'fuzzy_feature' property and simply add the
            # concerned feature to the features of the component.
            if not fuzzy_feature:
                comp.features.append(feat)
            else:
                comp.fuzzy_features.append(feat)


    def _handle_features_sbol_to_gb(self, seq_rec: SeqRecord, obj: sbol3.Component, seq_rec_features: List) -> None:
        """Helper function for resetting sequence features and their qualifiers to GenBank,
        by using a modified, extended SBOL3 Sequence Feature class - Feature_GenBank_Extension.
        :param seq_rec: GenBank SeqRecord instance for the record which contains sequnce features
        :param obj: Instance of the SBOL3 Component 
        :param seq_rec_features: List of features for the related sequence features
        """
        # parse if sbol object has any features
        if not obj.features:
            return
        feat_order = {}
        seq_rec_features = []
        # converting all sequence features
        for obj_feat in obj.features:
            # TODO: Also add ability to parse subcomponent feature type
            # Note: Currently we only parse sequence features from sbol3 to genbank
            if isinstance(obj_feat, sbol3.SequenceFeature):
                logging.info(
                    f"Parsing feature `{obj_feat.name}` for component `{obj.display_id}`"
                )
                # TODO: There may be multiple locations for a feature from sbol3; 
                #       add ability to parse them into a single genbank feature
                feat_loc_parts = []
                feat_loc_object = None
                feat_loc_positions = []
                feat_strand = self.BIO_STRAND_FORWARD
                for obj_feat_loc in obj_feat.locations:
                    feat_strand = self.BIO_STRAND_FORWARD
                    # feature strand value which denotes orientation of the location of the feature
                    # By default its 1 for SO_FORWARD orientation of sbol3 feature location, and -1 for SO_REVERSE
                    if obj_feat_loc.orientation == sbol3.SO_REVERSE:
                        feat_strand = self.BIO_STRAND_REVERSE
                    elif obj_feat_loc.orientation != sbol3.SO_FORWARD:
                        raise ValueError(f"Location orientation: `{obj_feat_loc.orientation}` for feature: \n \
                        `{obj_feat.name}` of component: `{obj.display_id}` is not a valid orientation.\n \
                        Valid orientations are `{sbol3.SO_FORWARD}`, `{sbol3.SO_REVERSE}`")
                    # TODO: Raise custom converter class ERROR for `else:`
                    # creating start and end Positions 
                    endPosition = ExactPosition(obj_feat_loc.end)
                    startPosition = ExactPosition(obj_feat_loc.start)
                    # if range, check for position being Before / After Positions
                    if isinstance(obj_feat_loc, sbol3.Range):
                        # change end and start Positions only if user has made integer entries into them
                        if obj_feat_loc.end_position != None:
                            endPosition = self.GENBANK_LOCATION_POSITION[obj_feat_loc.end_position](obj_feat_loc.end)
                        if obj_feat_loc.start_position != None:
                            startPosition = self.GENBANK_LOCATION_POSITION[obj_feat_loc.start_position](obj_feat_loc.start)
                    feat_loc_object = FeatureLocation(
                        start=startPosition,
                        end=endPosition,
                        strand=feat_strand,
                    )
                    feat_loc_parts.append(feat_loc_object)
                # sort feature locations lexicographically internally first
                # NOTE: If the feature location has an outer "complement" location operator, the sort needs to be in reverse order
                if obj_feat.orientation == sbol3.SO_REVERSE:
                    feat_loc_parts.sort(key=lambda loc: (loc.start, loc.end, loc.strand), reverse=True)
                else:
                    feat_loc_parts.sort(key=lambda loc: (loc.start, loc.end, loc.strand))
                for loc in feat_loc_parts:
                    feat_loc_positions += [loc.start, loc.end]
                if len(feat_loc_parts) > 1:
                    feat_loc_object = CompoundLocation(parts=feat_loc_parts, operator="join")
                elif len(feat_loc_parts) == 1:
                    feat_loc_object = feat_loc_parts[0]
                # action to perform if no location found?
                # else:

                # FIXME: order of features not same as original genbank doc?
                obj_feat_role = obj_feat.roles[0]
                # NOTE: The so2gb.csv data file has rows of format 'SO:xxxxxxx,<GenBank_Term>', 
                # and the obj_feat_role returns the URI (i.e 'https://identifiers.org/SO:xxxxxx').
                # The slicing and subtracting is done to obtain the 'SO:xxxxxxx' portion from the URI.
                obj_feat_role = obj_feat_role[
                    obj_feat_role.index(":", 6) - 2 :
                ]
                # Obtain sequence feature role from so2gb mappings
                feat_role = self.DEFAULT_GB_TERM
                if self.so2gb_map.get(obj_feat_role):
                    feat_role = self.so2gb_map[obj_feat_role]
                else:
                    logging.warning(f"Feature role: `{obj_feat_role}` for feature: `{obj_feat}` of component: \n \
                    `{obj.display_id}` has no corresponding ontology term for GenBank, using the default GenBank term, {self.DEFAULT_GB_TERM}")
                # create sequence feature object with label qualifier
                # TODO: create issue for presence of genbank file with features without the "label" qualifier
                # TODO: feat_strand value ambiguous in case of mulitple locations?
                feat = SeqFeature(
                    location=feat_loc_object, strand=feat_strand, type=feat_role
                )
                feat_order[feat] = feat_loc_positions
                if isinstance(obj_feat, self.Feature_GenBank_Extension):
                    keys = sorted(obj_feat.qualifier_key, key=lambda x: int(x.split(":", 1)[0]))
                    values = sorted(obj_feat.qualifier_value, key=lambda x: int(x.split(":", 1)[0]))
                    for qualifier_ind in range(len(keys)):
                        feat.qualifiers[keys[qualifier_ind].split(":", 1)[1]] = values[qualifier_ind].split(":", 1)[1]
                seq_rec_features.append(feat)

        # Sort features based on feature location start/end, lexicographically, and then by 
        # strand / number of qualifiers / type of feature string comparison
        seq_rec_features.sort(key=lambda feat: (feat_order[feat], feat.strand, len(feat.qualifiers), feat.type))
        seq_rec.features = seq_rec_features
