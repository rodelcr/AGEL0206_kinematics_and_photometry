-- Lua filter for pandoc -> latex:
-- Convert every inline Code element (markdown `...`) to a sequence
-- that wraps cleanly inside table cells. Strategy:
--   1. Escape the LaTeX special chars _ & # $ % ^ ~ (\, {, } also need
--      handling but rare in paths).
--   2. Insert \discretionary{}{}{} (an invisible break-opp) after
--      every / and after every escaped \_, so long paths break at
--      directory/word boundaries.
--   3. Wrap the result in \BreakablePath{...} which sets ttfamily.
--
-- This sidesteps the seqsplit/catcode problem: by the time our latex
-- sees the argument, special chars are properly escaped, so we don't
-- need to play catcode games.

local function escape_latex(s)
  -- Escape backslash first so it doesn't double-escape later
  s = s:gsub('\\', '\\textbackslash{}')
  s = s:gsub('([_&#%%$^~{}])', '\\%1')
  return s
end

function Code(el)
  if FORMAT ~= 'latex' and FORMAT ~= 'beamer' then
    return nil
  end
  local txt = escape_latex(el.text)
  -- Insert break opportunities after path/word separators
  txt = txt:gsub('(/)', '%1\\discretionary{}{}{}')
  txt = txt:gsub('(\\_)', '%1\\discretionary{}{}{}')
  txt = txt:gsub('(%-)', '%1\\discretionary{}{}{}')
  txt = txt:gsub('(%.)', '%1\\discretionary{}{}{}')
  return pandoc.RawInline('latex', '\\BreakablePath{' .. txt .. '}')
end
