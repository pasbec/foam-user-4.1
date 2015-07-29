#!/bin/bash
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Sourcing -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

PSD="$(dirname "${BASH_SOURCE[0]}")";
sourcePSD () { local psd="$PSD"; source "$1"; PSD="$psd"; }

# --------------------------------------------------------------------------- #

sourcePSD "$PSD/fct_error.sh"
sourcePSD "$PSD/fct_path.sh"
sourcePSD "$PSD/fct_dep.sh"

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

link_isLink ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return $?

    [ -L "$link" ] || return $?

    return 0
}


link_isNotLink ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return $?

    [ ! -L "$link" ] || return $?

    return 0
}


link_isLinkOrExit ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return $?

    if link_isLink "$link"; then
        return 0
    elif link_isNotLink "$link"; then
        error_fatal "File $link is no link." "$FUNCNAME" || return $?
    fi

    return -1
}


link_isNotLinkOrExit ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return $?

    if link_isNotLink "$link"; then
        return 0
    elif link_isLink "$link"; then
        error_fatal "File $link is no link." "$FUNCNAME" || return $?
    fi

    return -1
}


link_isDangling ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return $?

    link_isLinkOrExit "$link"|| return $?

    path_notExists "$(readlink $link)" || return $?

    return 0
}


link_isNotDangling ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return $?

    link_isLinkOrExit "$link"|| return $?

    path_exists "$(readlink $link)" || return $?

    return 0
}


link_rmDangling ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return $?

    if link_isDangling "$link"; then
        rm "$link" || return $?
    fi

    return 0
}


link_rmFindDangling ()
{
    local dir="$1"

    error_fatalIfEmptyVar 'dir' "$FUNCNAME" || return $?

    if dep_isHashedOrExit 'find'; then

        local linkList
            linkList="$(find "$dir" -type l)" || return $?

        local l
        for l in $linkList; do
            link_rmDangling "$l" || return $?
        done

    else
        return 1
    fi

    return 0
}


link_absPath ()
{
    local path="$1"

    error_fatalIfEmptyVar 'path' "$FUNCNAME" || return $?

    dep_isHashedOrExit 'python' && \
        python -c "import os.path; print os.path.abspath('$path')"

    [ $? -eq 0 ] || return $?

    return 0
}


link_relPath ()
{
    local path="$1"
    local refPath="${2:-.}"

    error_fatalIfEmptyVar 'path refPath' "$FUNCNAME" || return $?

    dep_isHashedOrExit 'python' && \
        python -c "import os.path; print os.path.relpath('$path','$refPath')"

    [ $? -eq 0 ] || return $?

    return 0
}


link_lnRel ()
{
    local dest="$1"
    local link="$2"

    error_fatalIfEmptyVar 'dest link' "$FUNCNAME" || return $?

    if path_exists "$dest"; then

        local linkDir
            linkDir="$(dirname "$link")" || return $?

        local destRelToLink
            destRelToLink="$(link_relPath "$dest" "$linkDir")" || return $?

        ln -s "$destRelToLink" "$link" || return $?

    else
        error_fatal "Link destination $dest not found." "$FUNCNAME" \
        || return $?
    fi

    return 0
}


link_lnAbsToRel ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return $?

    local linkDir
        linkDir="$(dirname $link)" || return $?

    local dest
        dest="$(readlink $link)" || return $?

    if path_isAbs "$dest"; then

        if path_exists "$dest"; then

            local destRelToLink
                destRelToLink="$(link_relPath $dest $linkDir)" || return $?

            if link_isNotLink "$link"; then

                error_fatal \
                    "File $link is no link." "$FUNCNAME" || return $?

            elif link_isLink "$link"; then

                echo "Relinking '$link' -> '$destRelToLink'"
                rm "$link" || return $?
                ln -s "$destRelToLink" "$link" || return $?

            fi

        elif path_notExists "$dest"; then

            error_fatal \
                "Link destination $dest not found." "$FUNCNAME" || return $?

        else
            return 1
        fi

    elif path_isNotAbs "$dest"; then

        error_info \
            "Skipping '$link' with relative destination '$dest'." \
            "$FUNCNAME" || return $?

        path_exists "$dest" || \
            error_fatal "Link destination $dest not found." \
            "$FUNCNAME" || return $?

    else
        return 1
    fi

    return 0
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
