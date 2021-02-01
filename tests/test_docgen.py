from pele_platform.Models import docgen


PATH = '../docs/build_docs/source/developers/yaml.rst'

def test_model_generator():
    render = docgen.render_yaml_parser()
    with open(PATH, "w") as f:
        f.write(render)
