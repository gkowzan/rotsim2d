if (($args.count -eq 1) -and ($args[0] -eq "force")) {
    rm requirements.txt dev-requirements.txt interactive-requirements.txt
}

python -m piptools compile --extra-index-url=http://127.0.0.1:4040 setup.cfg
python -m piptools compile dev-requirements.in
python -m piptools compile interactive-requirements.in
python -m piptools sync requirements.txt dev-requirements.txt interactive-requirements.txt
