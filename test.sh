function foo() {
    if false; then
        return 1;
    fi
}

foo &&
    echo check
