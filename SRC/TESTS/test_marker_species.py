"""Exercise the native namelist reader and pusher species override."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def main():
    binary = Path(sys.argv[1]).resolve()
    template = Path(sys.argv[2]).read_text()
    # Test omission separately from an explicit zero pair.
    template = re.sub(r'^\s*marker_(mass_g|charge_e)\s*=.*\n', '', template,
                      flags=re.MULTILINE)
    cases = [
        ('legacy', None, None, None),
        ('zero', '0', '0', None),
        ('AUG_D', '3.3436e-24', '1', None),
        ('ITER_D', '3.345228e-24', '1', None),
        ('ITER_T', '5.017842e-24', '1', None),
        ('ITER_e', '9.1091e-28', '-1', None),
        ('missing_charge', '3.345228e-24', None, 'Marker charge must be nonzero'),
        ('missing_mass', None, '1', 'Marker mass must be positive'),
        ('negative_mass', '-3.345228e-24', '1', 'Marker mass must be positive'),
        ('nan_mass', 'NaN', '1', 'Invalid marker mass'),
        ('infinite_charge', '3.345228e-24', 'Infinity', 'Invalid marker charge'),
    ]
    with tempfile.TemporaryDirectory(prefix='gorilla-species-') as work:
        root = Path(work)
        for name, mass, charge, error in cases:
            directory = root / name
            directory.mkdir()
            assignments = ''
            if mass is not None:
                assignments += f' marker_mass_g={mass},\n'
            if charge is not None:
                assignments += f' marker_charge_e={charge},\n'
            text = template.replace('&GORILLANML', '&GORILLANML\n' + assignments, 1)
            (directory / 'gorilla.inp').write_text(text.rstrip() + '\n')
            result = subprocess.run([str(binary), mass or '0', charge or '0'],
                                    cwd=directory, text=True, capture_output=True)
            output = result.stdout + result.stderr
            if error is None:
                assert result.returncode == 0, (name, output)
                assert 'PASS: marker species and pusher ratio' in output, (name, output)
            else:
                assert result.returncode != 0 and error in output, (name, output)
            print(f'PASS: {name}')


if __name__ == '__main__':
    main()
