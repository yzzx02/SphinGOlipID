import pytest

from scripts.build_sphingolipid_name_catalog import reference_type_name


@pytest.mark.parametrize('name,expected',[
    ('Cer(d18:1/16:0)','Cer'),
    ('M(IP)2C(t42:1)(2OH)','M(IP)2C'),
    ('1-O-acetyl-ω-linoleoyloxy-Cer(6-t100:1)','1-O-acetyl-ω-linoleoyloxy-Cer'),
])
def test_reference_type_is_separate_from_composition_and_modification_position(name,expected):
    assert reference_type_name(name)==expected
