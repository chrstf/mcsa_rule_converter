known_atom_attrs = {
    "id": "id",
    "x2": "x2",
    "y2": "y2",
    "elementType": "element_type",
    "mrvExtraLabel": "mrv_extra_label",
    "lonePair": "lone_pair",
    "formalCharge": "formal_charge",
    "mrvAlias": "mrv_alias",
    "hydrogenCount": "h_count",
    "sgroupRef": "sgroup_ref",
    "rgroupRef": "r_group_ref",
    "mrvMap": "mrv_map",
    "radical": "radical",
    "mrvPseudo": "mrv_pseudo",
    "mrvQueryProps": "mrv_query_props",
    "isotope": "isotope",
    "spinMultiplicity": "spin_multiplicity",
    "ligandOrder": "ligand_order",
    "attachmentOrder": "attachment_order",
    "attachmentPoint": "attachment_point",
    "sgroupAttachmentPoint": "sgroup_attachment_point"
}

known_bond_attrs = {
    "id": "id",
    "atomRefs2": "atom_refs_2",
    "atomRefs3": "atom_refs_3",
    "atomRefs4": "atom_refs_4",
    "order": "order",
    "convention": "convention",
    "bondStereo": "bond_stereo",
    "queryType": "query_type"
}

known_e_container_attrs = {
    "id": "id",
    "occupation": "occupation",
    "radical": "radical"
}

known_electron_attrs = {
    "atomRefs": "atom_refs",
    "difLoc": "dif_loc"
}

known_e_flow_attrs = {
    "arcAngle": "arc_angle",  # Angle of the electron flow arrow arc.
    "id": "id",
    "baseElectronContainerIndex": "b_e_container_index",  # Describes which electroncontainer holds electron
    "baseElectronIndexInContainer": "b_e_index_in_container",  # Describes which electron is base electron.
    "headFlags": "head_flags",
    "tailFlags": "tail_flags",
    "headWidth": "head_width",
    "headLength": "head_length",
    "headSkip": "head_skip",
    "tailSkip": "tail_skip"

}

known_atom_set_point_attrs = {
    "atomRefs": "atom_refs",
    "weights": "weights"
}

known_flow_base_point = {
    "atomRef": "atom_ref"
}