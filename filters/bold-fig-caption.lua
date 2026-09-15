-- In the "Figure legends" section, Quarto renders each caption as a Para
-- whose first inline is a raw openxml fragment setting the paragraph style
-- to "ImageCaption", followed by the resolved crossref label (e.g.
-- Str "Fig." Str "\160" Str "1" Str ":") and then the caption text itself.
-- This wraps just the label ("Fig. 1:") in Strong so it renders bold in the
-- docx output, without affecting inline citations like "(Fig. 1)" elsewhere
-- in the prose, which have no such raw style marker.
function Para(p)
  local content = p.content

  if
    #content == 0
    or content[1].t ~= "RawInline"
    or content[1].format ~= "openxml"
    or not content[1].text:find("ImageCaption")
  then
    return p
  end

  local label_end = nil
  for i = 2, math.min(#content, 6) do
    if content[i].t == "Str" and content[i].text == ":" then
      label_end = i
      break
    end
  end

  if not label_end then
    return p
  end

  local label_inlines = {}
  for i = 2, label_end do
    table.insert(label_inlines, content[i])
  end

  local new_content = { content[1], pandoc.Strong(label_inlines) }
  for i = label_end + 1, #content do
    table.insert(new_content, content[i])
  end

  p.content = new_content
  return p
end
