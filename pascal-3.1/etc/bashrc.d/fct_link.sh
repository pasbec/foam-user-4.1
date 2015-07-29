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

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return 1

    [ -L "$link" ] || return 1

    return 0
}


link_isLinkOrExit ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return 1

    link_isLink "$link" || error_fatal "File $link is no link." "$FUNCNAME"

    [ $? -eq 0 ] || return 1

    return 0
}


link_isNotLink ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return 1

    [ ! -L "$link" ] || return 1

    return 0
}


link_isNotLinkOrExit ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return 1

    link_isNotLink "$link" || error_fatal "File $link is no link." "$FUNCNAME"

    [ $? -eq 0 ] || return 1

    return 0
}


link_isDangling ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return 1

    link_isLinkOrExit "$link"|| return 1

    local dest
        dest=$(basename $link) || return 1
        dest=$(readlink $dest) || return 1

    path_notExists "$dest" || return 1

    return 0
}


link_isNotDangling ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return 1

    link_isLinkOrExit "$link"|| return 1

    local dest
        dest=$(basename $link) || return 1
        dest=$(readlink $dest) || return 1

    path_exists "$dest" || return 1

    return 0
}


link_rmDangling ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return 1

    link_isDangling "$link" && rm "$link"

    [ $? -eq 0 ] || return 1

    return 0
}


link_rmFindDangling ()
{
    local dir="$1"

    error_fatalIfEmptyVar 'dir' "$FUNCNAME" || return 1

    if dep_isHashedOrExit 'find'; then

        local linkList
            linkList="$(find "$dir" -type l)" || return 1

        local l
        for l in $linkList; do

            local link
                link="$l/$(basename $l)" || return 1

            link_rmDangling "$link" || return 1

        done

    else
        return 1
    fi

    return 0
}


link_absPath ()
{
    local path="$1"

    error_fatalIfEmptyVar 'path' "$FUNCNAME" || return 1

    dep_isHashedOrExit 'python' && \
        python -c "import os.path; print os.path.abspath('$path')"

    [ $? -eq 0 ] || return 1

    return 0
}


link_relPath ()
{
    local path="$1"
    local refPath="${2:-.}"

    error_fatalIfEmptyVar 'path refPath' "$FUNCNAME" || return 1

    dep_isHashedOrExit 'python' && \
        python -c "import os.path; print os.path.relpath('$path','$refPath')"

    [ $? -eq 0 ] || return 1

    return 0
}


link_lnRel ()
{
    local dest="$1"
    local link="$2"

    error_fatalIfEmptyVar 'dest link' "$FUNCNAME" || return 1

    path_exists "$dest" || \
        error_fatal "Link destination $dest not found." "$FUNCNAME"

    [ $? -eq 0 ] || return 1

    local linkDir
        linkDir="$(dirname "$link")" || return 1
    local destRelToLink
        destRelToLink="$(link_relPath "$dest" "$linkDir")" || return 1

    ln -s "$destRelToLink" "$link" || return 1

    return 0
}


link_lnAbsToRel ()
{
    local link="$1"

    error_fatalIfEmptyVar 'link' "$FUNCNAME" || return 1

    local linkDir
        linkDir="$(dirname $link)" || return 1
    local dest
        dest="$(readlink $link)" || return 1

    if path_isAbs "$dest"; then

        if path_exists "$dest"; then

            local destRelToLink="$(link_relPath $dest $linkDir)"

            if link_isNotLink "$link"; then

                error_fatal \
                    "File $link is no link." "$FUNCNAME" || return 1

            elif link_isLink "$link"; then

                echo "Relinking '$link' -> '$destRelToLink'"
                rm "$link" || return 1
                ln -s "$destRelToLink" "$link" || return 1

            fi

        elif path_notExists; then

            error_fatal \
                "Link destination $dest not found." "$FUNCNAME" || return 1

        else
            return 1
        fi

    elif path_isNotAbs; then

        error_info \
            "Skipping '$link' with relative destination '$dest'." \
            "$FUNCNAME" || return 1

        path_exists "$dest" || \
            error_fatal "Link destination $dest not found." \
            "$FUNCNAME" || return 1

    else
        return 1
    fi

    return 0
}

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
