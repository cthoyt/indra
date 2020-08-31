import logging
from collections import OrderedDict as _o


logger = logging.getLogger(__name__)

default_ns_order = ['WM', 'UN', 'HUME', 'SOFIA', 'CWMS']


class Concept(object):
    """A concept/entity of interest that is the argument of a Statement

    Parameters
    ----------
    name : str
        The name of the concept, possibly a canonicalized name.
    db_refs : dict
        Dictionary of database identifiers associated with this concept.
    """
    def __init__(self, name, db_refs=None):
        self.name = name
        self.db_refs = db_refs if db_refs else {}

    def matches(self, other):
        return self.matches_key() == other.matches_key()

    def matches_key(self):
        key = self.entity_matches_key()
        return str(key)

    def entity_matches(self, other):
        return self.entity_matches_key() == other.entity_matches_key()

    def entity_matches_key(self):
        # Get the grounding first
        db_ns, db_id = self.get_grounding()
        # If there's no grounding, just use the name as key
        if not db_ns and not db_id:
            return self.name
        return str((db_ns, db_id))

    def equals(self, other):
        matches = (self.name == other.name) and \
                  (self.db_refs == other.db_refs)
        return matches

    def get_grounding(self, ns_order=None):
        ns_order = ns_order if ns_order else default_ns_order
        for db_ns in ns_order:
            db_id = self.db_refs.get(db_ns)
            if not db_id:
                continue
            if isinstance(db_id, (list, tuple)):
                db_id = db_id[0]
                if isinstance(db_id, (list, tuple)):
                    db_id = sorted(db_id, key=lambda x: x[1],
                                   reverse=True)[0][0]
            return db_ns, db_id
        return None, None

    def isa(self, other, ontology):
        # Get the namespaces for the comparison
        (self_ns, self_id) = self.get_grounding()
        (other_ns, other_id) = other.get_grounding()
        # If one of the agents isn't grounded to a relevant namespace,
        # there can't be an isa relationship
        if not all((self_ns, self_id, other_ns, other_id)):
            return False
        # Check for isa relationship
        return ontology.isa(self_ns, self_id, other_ns, other_id)

    def is_opposite(self, other, ontology):
        # Get the namespaces for the comparison
        (self_ns, self_id) = self.get_grounding()
        (other_ns, other_id) = other.get_grounding()
        # If one of the agents isn't grounded to a relevant namespace,
        # there can't be an is_opposite relationship
        if not all((self_ns, self_id, other_ns, other_id)):
            return False
        # Check for is_opposite relationship
        return ontology.is_opposite(self_ns, self_id,
                                    other_ns, other_id)

    def refinement_of(self, other, ontology, entities_refined=False):
        # Make sure the Agent types match
        if type(self) != type(other):
            return False

        # Check that the basic entity of the agent either matches or is related
        # to the entity of the other agent. If not, no match.
        # If the entities, match, then we can continue
        if entities_refined:
            return True
        if self.entity_matches(other):
            return True
        if self.isa(other, ontology):
            return True
        return False

    def to_json(self):
        json_dict = _o({'name': self.name})
        json_dict['db_refs'] = self.db_refs
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        name = json_dict.get('name')
        db_refs = json_dict.get('db_refs', {})
        if not name:
            logger.error('Concept missing name.')
            return None
        # This fixes the fact that scored lists of groundings
        # are deserialized as lists of lists instead of lists
        # of tuples.
        for key, val in db_refs.items():
            if isinstance(val, list):
                db_refs[key] = [tuple(v) for v in val]
        concept = Concept(name, db_refs=db_refs)
        return concept

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)

