import os
import glob
import re
import pytest
from subprocess import Popen, PIPE


def execute_binary(binary_path, data_graph, query_graph):
    """Executes the orbitsi binary with the given graphs and graphlet size."""
    # The configuration under test comes from environment variables; unset ones leave the
    # binary's defaults.
    graphlet_size = os.environ.get('GRAPHLET_SIZE', '4')
    iterate_val = os.environ.get('ITERATE', None)
    iterate_flag = f' --iterate {iterate_val}' if iterate_val is not None else ''
    use_full_graph = os.environ.get('USE_FULL_GRAPH', '0')
    full_graph_flag = ' --use-full-graph' if use_full_graph == '1' else ''
    cbj = os.environ.get('CBJ', '0')
    cbj_flag = ' --cbj' if cbj == '1' else ''
    count_only = os.environ.get('COUNT_ONLY', '0')
    count_only_flag = ' --count-only' if count_only == '1' else ''
    arc_consistency = os.environ.get('ARC_CONSISTENCY', '1')   # the binary's default
    arc_consistency_flag = ' --no-arc-consistency' if arc_consistency == '0' else ''
    order_mode = os.environ.get('ORDER_MODE', None)   # unset = the binary's default
    order_mode_flag = f' --order-mode {order_mode}' if order_mode else ''

    extra_flags = os.environ.get('EXTRA_FLAGS', '')   # e.g. '--no-symmetry --no-edge-orbit-filter'
    extra_flag = f' {extra_flags}' if extra_flags else ''

    command = (f'{binary_path} --data {data_graph} --pattern {query_graph} '
               f'--graphlet-size {graphlet_size}{iterate_flag}{full_graph_flag}{cbj_flag}{count_only_flag}'
               f'{arc_consistency_flag}{order_mode_flag}{extra_flag}')
               
    process = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    std_output, std_error = process.communicate()
    return process.returncode, std_output.decode('utf-8'), std_error.decode('utf-8')


def parse_output(output):
    """Parses the output of the orbitsi binary to find the number of matches."""
    match = re.search(r"Matches found: (\d+)", output)
    return int(match.group(1)) if match else -1


@pytest.fixture(scope="session")
def project_root():
    """Get absolute path to project root (parent of 'test' directory)."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


@pytest.fixture(scope="session")
def binary_path(project_root):
    """Path to compiled orbitsi binary."""
    path = os.environ.get('ORBITSI_BIN') or os.path.join(project_root, "build", "orbitsi")
    if not os.path.isfile(path):
        pytest.skip(f"Binary not found at '{path}' — build it first.")
    return path


@pytest.fixture(scope="session")
def expected_results(project_root):
    """Load expected results from test/expected_output.res."""
    expected_file = os.path.join(project_root, "test", "expected_output.res")
    if not os.path.isfile(expected_file):
        pytest.skip(f"Expected results file not found at '{expected_file}'")

    results = {}
    with open(expected_file, 'r') as f:
        for line in f:
            if line.strip():
                name, count = line.strip().split(':')
                results[name.strip()] = int(count.strip())
    return results


@pytest.fixture(scope="session")
def data_graph_path(project_root):
    """Path to the main data graph."""
    path = os.path.join(project_root, "test", "data_graph", "HPRD.graph")
    if not os.path.isfile(path):
        pytest.skip(f"Data graph not found at '{path}'")
    return path


@pytest.fixture(scope="session")
def query_graphs(project_root):
    """List all query graph files."""
    folder = os.path.join(project_root, "test", "query_graph")
    paths = glob.glob(f'{folder}/*.graph')
    if not paths:
        pytest.skip(f"No query graphs found in '{folder}'")
    return paths


@pytest.mark.parametrize(
    "query_path",
    glob.glob(os.path.join(os.path.dirname(__file__), "query_graph", "*.graph")),
)
def test_query_correctness(binary_path, data_graph_path, query_path, expected_results):
    """Run the orbitsi binary for each query graph and verify match counts."""
    query_name = os.path.splitext(os.path.basename(query_path))[0]

    if query_name not in expected_results:
        pytest.skip(f"No expected result found for '{query_name}'")

    rc, std_out, std_err = execute_binary(binary_path, data_graph_path, query_path)

    assert rc == 0, f"Binary exited with non-zero status {rc}\nSTDERR:\n{std_err}"

    output_matches = parse_output(std_out)
    expected_matches = expected_results[query_name]

    assert output_matches == expected_matches, (
        f"Mismatch for '{query_name}': expected {expected_matches}, got {output_matches}"
    )